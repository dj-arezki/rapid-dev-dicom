# *****************************************************************
#
# IBM Confidential
# OCO Source Materials
# 
# (C) Copyright IBM Corp. 2019 All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
#
# The source code for this program is not published or otherwise
# divested of its trade secrets, irrespective of what has been
# deposited with the U.S. Copyright Office.
#
# *****************************************************************
# -*- coding: utf-8 -*-
"""
Class for reading and querying sets of DICOM images.

"""

import os
import numpy as np
import collections
from Dicom.DicomScanner import DicomScan
from Dicom.DicomFile import DicomFile
from Dicom.DicomImage import DicomImage


class DicomLoader(object):
    """
    Class for reading and querying DICOM series. See the __init__
    function for a description of the (python) attributes available
    after a series has been read.
    """

    def __init__(self, fileNames=None, slice_checker = None):
        """
        Initialize the object, and optionally read a series from file.

        Parameters
        ----------
        fileNames : str or sequence of strs
          The (optional) name of a single file to be read, name of directory from which
          files should be read, or list of file names to be read. See Load for details.
        """

        # Instance UIDs for each slice
        self.sopInstUIDs = None

        # Single (arbitrary) DicomImage from the series for examining attributes
        self.image = None

        # Image Position (Patient) from each slice
        self.imgPosPats = None

        # Rescale intercept and slope for each slice
        self.rescales = None

        # Padding values for each slice
        self.padVals = None

        # Contrast windows for each slice
        self.windows = None

        # Numpy array of raw pixel in native type without any rescaling
        self.pixels = None

        # Volume dimensions and spacing in z,y,x order
        self.nVoxels = None
        self.spacing = None
        self.origin  = None

        # 4x4 transformation matrix from grid space to DICOM patient space. Grid
        # space is the rectilinear space whose axes coincide with the grid axes,
        # whose origin is at the centre of the bottom front right voxel, and in
        # which voxels are separated by the voxel spacing
        self.gridToPat = None

        # Call Load if a filename specified
        if fileNames: self.Load(fileNames, slice_checker)

    # ========================================================================
    def Load(self, fileNames, slice_checker):
        """
        Load a DICOM series from file.

        Parameters
        ----------
        fileNames : str or sequence of strs
          The name of a single image file to be read, name of directory from which image
          files should be read, or list of image file names to be read. If a single file,
          the containing directory is scanned for image files with the same series UID,
          and they are all loaded. If a directory, it is scanned for image files, and the
          set of files belonging to the largest series is loaded. If a list of image file
          names with the same series UID, they are all loaded.
        """

        # Make sure the input is a path or a list of paths
        if not isinstance(fileNames, str) and (not isinstance(fileNames, collections.
            Sequence) or len([p for p in fileNames if not isinstance(p, str)]) != 0):
            raise TypeError("Expected a path or a sequence of paths")

        #-----------------------
        # Identify files to load
        #-----------------------

        # If it's a single path
        if isinstance(fileNames, str):

            # Initialize
            ser2LoadUid = None
            serTag = (0x0020,0x000E)

            # If the path is to a file
            if os.path.isfile(fileNames):

                # Load the file and get the series UID
                dcmFile = DicomFile(fileNames)
                ser2LoadUid = dcmFile.GetAttribute(serTag)

                # Reset to point at containing directory
                fileNames = os.path.dirname(fileNames)

            # Get the list of files in the directory
            dir = os.path.dirname(fileNames)
            pathNames = []
            for root, dir, file_list in os.walk(dir):
                for file in file_list:
                    pathNames.append(os.path.join(root, file))

            # Collect all paths and series UIDs from files
            ser2Path = {}
            def SerPathCollector(path, tvDict):
                if serTag in tvDict:
                    serUid = tvDict[serTag]
                    if serUid in ser2Path:
                        ser2Path[serUid].append(path)
                    else: ser2Path[serUid] = [path]
            DicomScan(pathNames, (serTag,), SerPathCollector)

            # Make sure we found something
            if not ser2Path: raise IOError("No DICOM files found")

            # If we don't have a series UID yet
            if not ser2LoadUid:

                # Pick the one with the most files
                maxLen = 0
                for uid in ser2Path:
                    if maxLen < len(ser2Path[uid]):
                        maxLen = len(ser2Path[uid])
                        ser2LoadUid = uid

            # Get list of files to load
            fileNames = ser2Path[ser2LoadUid]

        #------------------------------------------
        # Check attribute consistency across slices
        #------------------------------------------

        # Attributes we expect to be the same across images
        sameImgTags = (
            (0x0008,0x0060), #  0 - Modality
            (0x0008,0x0008), #  1 - Image Type
            (0x0028,0x0010), #  2 - Rows
            (0x0028,0x0011), #  3 - Columns
            (0x0028,0x0002), #  4 - Samples Per Pixel
            (0x0028,0x0100), #  5 - Bits Allocated
            (0x0028,0x0101), #  6 - Bits Stored
            (0x0028,0x0102), #  7 - High Bit
            (0x0028,0x0103), #  8 - Pixel Representation
            (0x0028,0x0030), #  9 - Pixel Spacing
            (0x0028,0x0004), # 10 - Photometric Interpretation
            (0x0020,0x000D), # 11 - Study Instance UID
            (0x0020,0x000E), # 12 - Series Instance UID
            (0x0020,0x0052), # 13 - Frame Of Reference UID
            (0x0018,0x0050), # 14 - Slice Thickness
        )

        # Dictionary to compare across slices, and flag
        # recording whether dictionaries are all the same
        imgDict = {}
        notConsistent = []

        # Check consistency with scan
        def CheckConsistency(path, tvDict):
            if not imgDict: imgDict.update(tvDict)
            elif imgDict != tvDict: notConsistent.append(True)
        DicomScan(fileNames, sameImgTags, CheckConsistency)
        if notConsistent: raise RuntimeError("Images are inconsistent")

        # Error checking of these attributes
        if sameImgTags[2] not in imgDict or imgDict[sameImgTags[2]] <= 0:
            raise RuntimeError("Rows absent or non-positive")

        if sameImgTags[3] not in imgDict or imgDict[sameImgTags[3]] <= 0:
            raise RuntimeError("Columns absent or non-positive")

        if sameImgTags[4] not in imgDict or imgDict[sameImgTags[4]] != 1:
            raise RuntimeError("Samples Per Pixel absent or not equal to 1")

        if sameImgTags[5] not in imgDict or imgDict[sameImgTags[5]] not in (8, 16):
            raise RuntimeError("Bits Allocated absent or not 8 or 16")

        if sameImgTags[6] not in imgDict or imgDict[sameImgTags[6]] <= 0:
            raise RuntimeError("Bits Stored absent or non-positive")

        if sameImgTags[7] not in imgDict or imgDict[sameImgTags[7]] <= 0:
            raise RuntimeError("High Bit absent or non-positive")

        if sameImgTags[8] not in imgDict or imgDict[sameImgTags[8]] not in (0, 1):
            raise RuntimeError("Pixel Representation absent or not 0 or 1")

        if (sameImgTags[9] not in imgDict or len(imgDict[sameImgTags[9]]) != 2 or
            imgDict[sameImgTags[9]][0] <= 0.0 or imgDict[sameImgTags[9]][1] <= 0.0):
            raise RuntimeError("Pixel Spacing absent or non-positive")

        if sameImgTags[10] not in imgDict or imgDict[sameImgTags[10]] != "MONOCHROME2":
            raise RuntimeError("Photometric Interpretation absent or not MONOCHROME2")

        #------------------------------------------
        # Load attributes that differ across slices
        #------------------------------------------

        # Initialize persistent per-slice collections
        self.sopInstUIDs = []
        self.image = None
        self.imgPosPats = []
        self.rescales = []
        self.padVals = []
        self.windows = []
        self.pixels = None

        # To make sure we don't load the same slice twice
        instUidMap = {}

        # For image orientation
        imgOrientPats = []

        # For each slice of data
        pixels = []

        # Loop over slices
        for fn in fileNames:

            # Try to load the image
            img = DicomImage(fn)

            # Get the instance UID
            instUid = img.GetAttribute((0x0008,0x0018))
            if not instUid: raise RuntimeError("SOP Instance UID absent")

            # Make sure we haven't loaded it already
            if instUid in instUidMap: continue
            else: instUidMap[instUid] = None

            # Save the UID
            self.sopInstUIDs.append(instUid)

            # Save the image if not done yet
            if not self.image: self.image = img

            # Get the image position (patient)
            imgPosPat = img.GetAttribute((0x0020,0x0032))
            if not imgPosPat: raise RuntimeError("Image Position (Patient) absent")
            if len(imgPosPat) != 3: raise RuntimeError("Image Position (Patient) invalid")
            self.imgPosPats.append([float(ipp) for ipp in imgPosPat])

            # Get the image orientation (patient)
            imgOrientPat = img.GetAttribute((0x0020,0x0037))
            if not imgOrientPat: raise RuntimeError("Image Orientation (Patient) absent")
            if len(imgOrientPat) != 6: raise RuntimeError("Image Orientation (Patient) invalid")
            imgOrientPats.append([float(iop) for iop in imgOrientPat])

            # Get the image rescale intercept and slope
            rescale = [img.GetAttribute((0x0028,0x1052)),img.GetAttribute((0x0028,0x1053))]

            # Only mandatory for CT and PET
            modality = img.GetAttribute((0x0008,0x0060))
            if modality == "CT":

                # Save the mandatory values
                self.rescales.append(
                    [float(rescale[0]) if rescale[0] else 0.0,
                     float(rescale[1]) if rescale[1] else 1.0])

            elif modality == "PT":

                # Save the mandatory values
                if not rescale[0] or not rescale[1]: raise \
                    RuntimeError("Rescale Intercept or Slope absent")
                self.rescales.append([float(rescale[0]), float(rescale[1])])

            else:

                # Save any optional values
                self.rescales.append([0.0,1.0])
                if rescale[0]: self.rescales[-1][0] = float(rescale[0])
                if rescale[1]: self.rescales[-1][1] = float(rescale[1])

            # Get the pixel padding value, if any
            padVal = img.GetAttribute((0x0028,0x0120))
            self.padVals.append(padVal)

            # Get the contrast windows, if any
            windCtr = img.GetAttribute((0x0028,0x1050))
            windWid = img.GetAttribute((0x0028,0x1051))
            if isinstance(windCtr, str) and isinstance(windWid, str):

                self.windows.append([[float(windCtr),float(windWid)]])

            elif (isinstance(windCtr, collections.Sequence) and
                  isinstance(windWid, collections.Sequence) and
                  len(windCtr) == len(windWid)):

                self.windows.append([])
                for i in range(len(windCtr)):
                    self.windows[-1].append([float(windCtr[i]), float(windWid[i])])

            else: self.windows.append([])

            if slice_checker is not None:
                slice_checker(img)

            # Save the pixel data
            pixels.append(img.GetPixelData(native=True))

        #-----------------------------------------
        # Sort images and validate volume geometry
        #-----------------------------------------

        # The maximum amount by which unit vectors can deviate
        # from having unit norm
        unitVecMagDev = 0.02

        # The difference between the minimum and maximum gap
        # size must not be greater than this amount
        maxAllowedGapDev = 0.02; # mm

        # The direction cosine between vectors separating adjacent pairs
        # of slices and between these vectors and the acquisition axis must
        # not differ from 1 by more than this amount
        maxAllowedDirCosDev = 0.0075; # 7.0 degrees

        # The fractional difference between the lengths of
        # vectors separating slices must not exceed this amount
        maxAllowedGapFractDev = 0.05; # 1 part in 20

        # Get the orientation and acquisition vectors
        imgOrientPats = np.array(imgOrientPats, dtype=np.float32)
        acqVecs = np.cross(imgOrientPats[:,:3], imgOrientPats[:,3:])

        # Validate orientation data
        if any(np.fabs(np.linalg.norm(imgOrientPats[:,:3], axis=1) - 1.0) > unitVecMagDev):
            raise RuntimeError("Image Orientation (Patient) row vector norm not 1")

        if any(np.fabs(np.linalg.norm(imgOrientPats[:,3:], axis=1) - 1.0) > unitVecMagDev):
            raise RuntimeError("Image Orientation (Patient) column vector norm not 1")

        if any(np.fabs(np.linalg.norm(acqVecs, axis=1) - 1.0) > unitVecMagDev):
            raise RuntimeError("Image Orientation (Patient) row and column vectors not orthogonal")

        if any(np.fabs(np.dot(acqVecs[0,:], acqVecs.T) - 1.0) > maxAllowedDirCosDev):
            raise RuntimeError("Image Orientation (Patient) indicates slices not parallel")

        # Get the distances along the acquisition vector from the first image
        imgPosPats = np.array(self.imgPosPats, dtype=np.float32)
        imgPosPats -= imgPosPats[0,:] * np.ones((imgPosPats.shape[0],1), dtype=np.float32)
        imgPosPats = np.dot(imgPosPats, acqVecs[0,:].T)

        # Get the re-ordering permutation
        ordPerm = np.argsort(imgPosPats)

        if len(imgPosPats) <= 1:
            raise RuntimeError("Too few slices in series")

        # Validate the gap spacings
        imgPosPats = imgPosPats[ordPerm]
        gaps = np.diff(imgPosPats)
        maxGap = gaps.max()
        minGap = gaps.min()

        if minGap < 1.0e-6:
            raise RuntimeError("One or more slices are co-incident")

        if maxGap - minGap > maxAllowedGapDev:
            raise RuntimeError("Inhomogeneous gap spacing between "
                "slices: %g to %g" % (gaps.min(), gaps.max()))

        # Re-order the per-slice arrays
        self.sopInstUIDs = [self.sopInstUIDs[i] for i in ordPerm]
        self.imgPosPats = [self.imgPosPats[i] for i in ordPerm]
        self.rescales = [self.rescales[i] for i in ordPerm]
        self.padVals = [self.padVals[i] for i in ordPerm]
        self.windows = [self.windows[i] for i in ordPerm]
        pixels = [pixels[i] for i in ordPerm]

        # Initialize pixel data
        nX = imgDict[sameImgTags[3]]
        nY = imgDict[sameImgTags[2]]
        nZ = len(pixels)
        self.pixels = np.empty((nZ,nY,nX), dtype=pixels[0].dtype)

        # Loop over each slice and transfer slice pixels to volume
        for i in range(len(pixels)): self.pixels[i,:,:] = pixels[i][:,:]

        # Set the volume dimensions and spacing in x,y,z order
        dX = imgDict[sameImgTags[9]][1] # Column spacing
        dY = imgDict[sameImgTags[9]][0] # Row spacing
        dZ = np.median(gaps) # Median slice gap spacing

        self.nVoxels = self.pixels.shape
        self.spacing = [dZ, dY, dX]
        self.origin = self.imgPosPats[0][::-1]

        # Set the transformation from grid space to patient space
        self.gridToPat = np.identity(4, dtype=np.float32)
        self.gridToPat[0:3,0] = imgOrientPats[0,:3]
        self.gridToPat[0:3,1] = imgOrientPats[0,3:]
        self.gridToPat[0:3,2] = acqVecs[0,:]
        self.gridToPat[0:3,3] = self.imgPosPats[0]

    def GetImage(self):
        """
        Return a representative image from the acquisition.
        """
        return self.image

    def GetDimensions(self):
        """
        Return the dimensions of the acquisition.
        """
        return self.nVoxels

    def GetSpacing(self):
        """
        Return the spacing of the acquisition.
        """
        return self.spacing

    def GetOrigin(self):
        """
        Return the origin of the acquisition.
        """
        return self.origin

    def GetImagePositionPatients(self):
        """
        Return the Image Position (Patient) attributes
        from the acquisition in [Z,Y,X] order.
        """
        return [ipp[::-1] for ipp in self.imgPosPats]

    def GetSOPInstanceUIDs(self):
        """
        Return the SOP Instance UIDs from the acquisition.
        """
        return self.sopInstUIDs

    def GetRescaleInterceptAndSlopes(self):
        """
        Return the Rescale Intercept and Rescale Slope attributes from the acquisition.
        """
        return self.rescales

    def GetPixelPaddingValues(self):
        """
        Return the list of pixel padding values from the acquisition.
        """
        return self.padVals

    def GetContrastWindows(self):
        """
        Return the set of contrast windows from the acquisition.
        """
        self.windows

    def GetPixelData(self, native = False):
        """
        Get a numpy array with the series pixel data.

        Parameters
        ----------
        native : bool
          Should the pixel data be returned in its native format? If False
          (the default), the pixel data will be returned as numpy.float32,
          after applying the rescale slope and intercept if they are
          present.
        """

        # Return to caller as appropriate type
        if not native:

            # Get rescaling intercept and slope arrays
            icept = np.array([s[0] for s in self.rescales], dtype=np.float32)
            slope = np.array([s[1] for s in self.rescales], dtype=np.float32)

            # Reshape for broadcasting
            icept = icept[:,np.newaxis,np.newaxis]
            slope = slope[:,np.newaxis,np.newaxis]

            # Return the rescaled array
            return slope * self.pixels.astype(np.float32) + icept

        else: return self.pixels
        
    def IsOrientedToPatient(self):
        """
        Returns true if the data is already oriented to the DICOM
        standard patient-based coordinate system, otherwise false.
        """
        return bool(self.gridToPat[0,0] > 0.9
                and self.gridToPat[1,1] > 0.9
                and self.gridToPat[2,2] > 0.9)

    def ReorientToPatient(self, pixels):
        """
        Re-orient an image volume to the DICOM standard patient-based 
        coordinate system, with X from right to left, Y from anterior 
        to posterior, and Z from inferior to superior. Raises RuntimeError
        exception on failure.

        Parameters
        ----------
        pixels : 3D numpy array
            Volume to be re-oriented according to the orientation 
            data associated with the instance's loaded DICOM volume.

        Returns
        -------
        pixels : 3D numpy array
            The re-oriented volume.
        spc : list with three floats
            The re-oriented spacing.
        org : list with three floats
            The re-oriented origin.
        origToPat : int32 numpy array with shape (3,4)
            Matrix for transforming pixel coordinates from
            the original space to the patient space.
        patToOrig : int32 numpy array with shape (3,4)
            Matrix for transforming pixel coordinates from
            patient space back to the original space.

        Usage example for origToPat and patToOrig
        -----------------------------------------
        Let xo, yo, zo be the coordinates of a voxel in the original
        volume, and let xp, yp, zp be the coordinates of the same voxel
        in the volume re-oriented to DICOM patient space. To convert
        one to the other:

            vOrig = np.array([zo, yo, xo], np.int32)
            vPat = origToPat[:,:3] @ vOrig + origToPat[:,3]
            zp, yp, xp = vPat.tolist()

            vPat = np.array([zp, yp, xp], np.int32)
            vOrig = patToOrig[:,:3] @ vPat + patToOrig[:,3]
            zo, yo, xo = vOrig.tolist()
        """

        # Do a quick initial check to make sure that the 
        # transformation will be a permutation of the axes.
        # If not, then raise exception
        dirCosThresh = 0.9
        if any(sum != 1 for sum in np.sum(np.where(np.fabs(
            self.gridToPat[:3,:3]) > dirCosThresh, 1, 0), -1)):
            raise RuntimeError("Patient orientation invalid")

        # This section determines the dominant patient space axis and its direction for 
        # each of the initial image space axes, i.e. which of the image row, column, and 
        # acquisition axes lie most closely along the patient x, y, and z axes. Knowing
        # which patient axis is initially closest to each image axis helps determine how
        # the image volume must be rotated so that the desired output is achieved, in 
        # which the row axis is along the patient positive x axis, the column axis is 
        # along the positive y axis, and the acquisition axis is along the positive z
        # axis. 
        # 
        # For example, if the current row axis is close to the patient positive z axis, 
        # then one constraint on the transformation that is implemented is that it should 
        # map the positive z axis to the positive x axis. The transformation is determined 
        # by establishing this mapping for each image axis

        # Copy and rename the gridToPat transform gridToImg since it is the transformation 
        # from the image grid to the image space, which may or may not be the same as patient 
        # space - that is to be determined below
        gridToImg = self.gridToPat

        # Array of direction cosine magnitudes and the image/patient axis pair they
        # describe, sorted into decreasing order of direction cosine magnitude
        dirCosImgPat = sorted([(np.fabs(gridToImg[j,i]), i, j) for i in 
            range(3) for j in range(3)], key=lambda k: k[0], reverse=True)

        # Array indicating which image axes have already been assigned
        # to a patient axis, and the patient axis to which they were 
        # assigned (X,Y,Z)
        imgPatAxis = [None]*3

        # Array holding the sign of the assignment 
        # from image axis to patient axis
        imgPatSign = [1]*3

        # Flags indicating which patient axes X,Y,Z 
        # have already been claimed by some image axis
        patAxisTaken = [False]*3

        # Assignment loop. Image axes are assigned to patient axes in order 
        # of decreasing direction cosine magnitude from the sort above. Loop 
        # executes if the image axis i hasn't been assigned a patient axis yet, 
        # and the patient axis dirCosImgPat[i].get<PA>() is available, dc is
        # the dir cosine, ia the image axis, pa the patient axis
        for dc, ia, pa in dirCosImgPat:
            if imgPatAxis[ia] is None and not patAxisTaken[pa]:

                # Assign the patient axis to the 
                # image axis and mark it taken
                imgPatAxis[ia] = pa
                patAxisTaken[pa] = True

                # Record the sign of the mapping
                if gridToImg[pa,ia] < 0.0:
                    imgPatSign[ia] = -1

        # Return failure if any axis unassigned
        if any(ipa is None for ipa in imgPatAxis):
            raise RuntimeError("Patient orientation invalid")

        # Initialize forward/backward pixel coordinate
        # transformations that will be returned
        origToPat = np.zeros((3,4), np.int32)
        patToOrig = np.zeros((3,4), np.int32)

        # Return default if no re-orientation is needed
        if all(ipa == i and ips == 1 for i, (ipa, ips) 
            in enumerate(zip(imgPatAxis, imgPatSign))):
            for i in range(3):
                origToPat[i,i] = 1
                patToOrig[i,i] = 1
            return pixels, self.spacing.copy(), \
                self.origin.copy(), origToPat, \
                    patToOrig

		# Construct the transformation 
        # matrix from image to patient space
        imgPatXform = np.zeros((3,3), np.int32)
        for i, (ipa, ips) in enumerate(zip(imgPatAxis, imgPatSign)):
            imgPatXform[ipa, i] = ips

        # Image dimensions and spacing and
        # patient dimensions and spacing
        nImg = self.nVoxels[::-1]
        dImg = self.spacing[::-1]
        nPat = [0]*3; dPat = [0.0]*3

		# Compute the patient space dimensions nPat, 
		# spacing dPat, and patient axes gridToPat
        gridToPat = np.zeros((3,3), np.float32)
        for i in range(3):
            for j in range(3):
                nPat[i] += imgPatXform[i,j] * nImg[j]
                dPat[i] += imgPatXform[i,j] * dImg[j]
                for k in range(3):
                    gridToPat[i,j] += gridToImg[i,k] * imgPatXform[j,k]

		# Take the absolute values of the dimensions and spacing, and 
		# compute offsets imgPatOffset for dimensions that were negative 
		# to keep them in the original range
        imgPatOffset = [0]*3
        for i in range(3):
            if nPat[i] < 0:
                nPat[i] = -nPat[i]
                dPat[i] = -dPat[i]
                imgPatOffset[i] = nPat[i] - 1

        # Compute the delta from the origin in image space to the origin 
        # in patient space, by computing the coordinates in image space 
        # of the (0,0,0) coordinate origin in patient space
        imgOrigToPatOrigDelta = [0]*3
        for i in range(3):
            for j in range(3):
                imgOrigToPatOrigDelta[i] -= \
                    imgPatXform[j,i] * imgPatOffset[j]

        # Use the position of the (0,0,0) coordinate origin in image space and 
        # the coordinate delta from the image space origin to the patient space 
        # origin to compute the position of the patient (0,0,0) coordinates in 
        # patient space
        oPat = gridToImg[:3,3].tolist()
        for i in range(3):
            for j in range(3):
                oPat[i] += gridToImg[i,j] * \
                    dImg[j] * imgOrigToPatOrigDelta[j]

        # Set forward and reverse transforms
        for i in range(3):
            for j in range(3):
                origToPat[i,j] = imgPatXform[2-i,2-j]
                patToOrig[i,j] = imgPatXform[2-j,2-i]
            origToPat[i,3] = imgPatOffset[2-i]
        for i in range(3):
            for j in range(3):
                patToOrig[i,3] -= patToOrig[i,j] * origToPat[j,3]

        # Perform any necessary axis flips
        if imgPatSign[2] == -1: pixels = pixels[::-1,::,::]
        if imgPatSign[1] == -1: pixels = pixels[::,::-1,::]
        if imgPatSign[0] == -1: pixels = pixels[::,::,::-1]

        # Perform any necessary permutations
        if any(ipa != i for i, ipa in enumerate(imgPatAxis)):
            if imgPatAxis[2] == 2:
                pixels = np.transpose(pixels, (1,0,2))
            elif imgPatAxis[2] == 1:
                if imgPatAxis[1] == 2:
                    pixels = np.transpose(pixels, (1,0,2))
                else:
                    pixels = np.transpose(pixels, (0,2,1))
                    pixels = np.transpose(pixels, (1,0,2))
            else:
                if imgPatAxis[1] == 2:
                    pixels = np.transpose(pixels, (1,0,2))
                    pixels = np.transpose(pixels, (0,2,1))
                else:
                    pixels = np.transpose(pixels, (2,1,0))

        # Make sure data remain contiguous and aligned
        pixels = np.require(pixels, requirements=["C","A"])

        # Return results in Z,Y,X order along with transformations
        return pixels, dPat[::-1], oPat[::-1], origToPat, patToOrig
