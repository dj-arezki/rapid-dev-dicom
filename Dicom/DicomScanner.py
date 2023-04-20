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
Function and class for gathering attributes from DICOM files.
"""

import collections
import os
import sys
import weakref as wr
from pathlib import Path

import _pickle as pickle
import gdcm
import numpy as np

from Dicom.DicomBase import (DicomAcquisition, DicomBase, DicomPatient,
                             DicomSeries, DicomStudy)
from Dicom.DicomConstants import *


def DicomScan(paths, groupElem, func):
    """
    Scan and gather attributes from DICOM files.

    Parameters
    ----------
    paths : str or sequence of str
        Names of files and directories to be scanned. Can be a
        mixture of both.
    groupElem: sequence of (int,int) tuples
        DICOM attributes to be extracted from all scanned files.
    func : function taking two arguments
        This function is called once for each file processed. The
        first argument is the path to the file, while the second is
        a dict mapping (group,elem) DICOM attribute tags to their
        associated value.
    """

    # Make sure the input is a path or a list of paths
    if not isinstance(paths, str) and (not isinstance(paths, collections.
        Sequence) or len([p for p in paths if not isinstance(p, str)]) != 0):
        raise TypeError("Expected a path or a sequence of paths")

    # Turn single path into list
    if isinstance(paths, str): paths = [paths]

    # Handle edge cases gracefully
    if len(paths) == 0: return None

    # Is the argument a group/element tuple?
    def isGroupElem(ge): return (isinstance(ge, tuple) and len(ge) == 2
        and isinstance(ge[0], int) and isinstance(ge[1], int))

    # Make sure the input is a group/elem or a list of group/elem tuples
    if not isGroupElem(groupElem) and (not isinstance(groupElem, collections.
        Sequence) or len([ge for ge in groupElem if not isGroupElem(ge)]) != 0):
        raise TypeError("Expected a list of group/elem tuples")

    # Turn single path into list
    if isGroupElem(groupElem): groupElem = [groupElem]

    # Handle edge cases gracefully
    if len(groupElem) == 0: return None

    # List of files we are going to look at
    files = []

    # Extensions we are not interested in
    noExt = ["txt","csv","bin","xml","zip","nii","gz","htm","html","exe","c","cpp","h","hdf5",
        "py","pyc","doc","docx","xls","xlsx","xlsm","ppt","pptx","jpg","bmp","img","png"]

    # Mapping from VR to Python type function
    vrMap = { "US":int, "UL":int, "SS":int, "SL":int,
        "FL":float, "FD":float, "DS":float, "IS":int }

    # Construct list of files
    for path in paths:

        # Skip if it doesn't exist
        if not os.path.exists(path): continue

        # If it's a directory
        if os.path.isdir(path):
            for root, dir, file_list in os.walk(path):
                for file in file_list:
                    path_in_dir = os.path.join(root, file)
                    if Path(path_in_dir).suffix not in noExt:
                        files.append(path_in_dir)

        elif os.path.isfile(path):
            # Add the single file
            files.append(path)

    # Turn off warnings
    gdcm.Trace.WarningOff()

    # Read dictionary to get VRs
    globInst = gdcm.Global.GetInstance()
    globInst.LoadResourcesFiles()
    dicts = globInst.GetDicts()

    # Create scanner
    scanPtr = gdcm.Scanner.New();
    scan = scanPtr.__ref__()

    # Add the attributes we are interested in
    for ge in groupElem: scan.AddTag(gdcm.Tag(ge[0], ge[1]));

    # Execute the scan
    success = scan.Scan(files)

    # Abort if we failed
    if not success: return

    # Check results for each file
    for fnam in files:

        # Declare empty dict
        tagValueDict = {}

        # Start iteration over the attributes
        pttv = gdcm.PythonTagToValue(scan.GetMapping(fnam))

        # Iterate until the end
        pttv.Start()
        while not pttv.IsAtEnd():

            # Get current tag and its VR
            tag = pttv.GetCurrentTag()
            vr = str(dicts.GetDictEntry(tag).GetVR())

            # Initialize value
            value = None

            # If not a sequence
            if vr != "SQ":

                # Get the value
                value = pttv.GetCurrentValue().strip()

                # Parse attribute contents
                if value:

                    # Split into component parts
                    value = value.split("\\")

                    # Separate and parse to appropriate type
                    if vr in vrMap:

                        try:

                            # Necessary to wrap this in try/except block
                            # as very rarely, the VR rules are violated
                            value = [vrMap[vr](attr) for attr in value]

                        except Exception:

                            # Advance and skip the rest, as
                            # if the tag was not present
                            pttv.Next()
                            continue

                    value = value if len(value) > 1 else value[0]

            else: value = "SQ" # Nested attributes not supported yet

            # Store in dictionary with tag
            tagValueDict[(tag.GetGroup(),tag.GetElement())] = value

            # Advance to next tag/value
            pttv.Next()

        # Call user-supplied function
        func(fnam, tagValueDict)

class DicomScanner(object):
    """
    Scan and gather attributes from DICOM files in a patient/study
    /series/acquisition hierarchy. See the classes in DicomBase for
    details of the objects generated by the scan.

    For the purposes of this class, an acquisition is a group of
    objects having the same Series Instance UID (0020,000E), and
    the same values of the following attributes:

        Rows                        (0028,0010)
        Columns                     (0028,0011)
        Pixel Spacing               (0028,0030)
        ImageType                   (0008,0008)
        Protocol Name               (0018,1030)
        Image Orientation (Patient) (0020,0037)

        Note that this does *not* include the acquisition number.

    NOTE: The output is NOT a deep copy! If you change it, you may
    corrupt this object's internal data structures! Consider it
    read-only.

    A typical loop to traverse all of the collected data might look
    like this:

    # Create the scanner and execute the scan
    scanner = DicomScanner()
    scanner.Scan("myPath")

    # Loop over all patients
    for pat in scanner:

        # Loop over all studies in patient
        for std in pat:

            # Loop over all series in study
            for ser in std:

                # Only interested in CT
                if ser.GetModality() != "CT": continue

                # Loop over all acquisitions in series
                for acq in ser:

                    # Get the list files for this acquisition
                    acqFiles = acq.GetFiles()

                    # Load them into a DicomSeries
                    dicomSer = DicomLoader(acqFiles)

                    # Get the pixel data as a numpy array
                    volume = dicomSer.GetPixelData()

                    ...

    """

    def __init__(self, dscFilePath=None):
        """
        Initialize the instance with attributes to be collected at each level of the
        hierarchy. Additional attributes to be queried may be provided at each level
        by calling the Add???Attribute() functions. See the __init__ function for the
        list of attributes at each level that are collected by default.
        """

        # Set the patient level tags
        self.patientTags = [
            PatientsNameTag,
            PatientsAgeTag,
            PatientsSexTag,
            PatientCommentsTag]

        # Set the study level tags
        self.studyTags = [
            StudyDescriptionTag,
            StudyDateTag,
            AccessionNumberTag]

        # Set the series level tags
        self.seriesTags = [
            SeriesDescriptionTag,
            SeriesModalityTag,
            ManufacturerTag,
            ManufacturerModelTag,
            FrameOfReferenceUIDTag,
            SOPClassUIDTag]

        # Set the acquisition level tags
        self.acquisitionTags = [
            ImageTypeTag,
            ImageCommentsTag,
            ProtocolNameTag,
            AcquisitionDateTag,
            AcquisitionTimeTag]

        # If we got a scanner file path
        if dscFilePath is not None:

            # Load it
            self.Load(dscFilePath)

        else: # No path provided

            # Initialize patients
            self.patients = []

    def __iter__(self):
        """
        Provide direct iterable access to list of patients.
        """
        return iter(self.patients)

    def AddPatientAttribute(self, groupElem):
        """
        Add the specified attribute to the set of patient attributes to be scanned.

        Parameters
        ----------
        groupElem : (int,int) tuple
            (group,elem) tag for the attribute to be added.
        """
        self.patientTags.append(groupElem)

    def AddStudyAttribute(self, groupElem):
        """
        Add the specified attribute to the set of study attributes to be scanned.

        Parameters
        ----------
        groupElem : (int,int) tuple
            (group,elem) tag for the attribute to be added.
        """
        self.studyTags.append(groupElem)

    def AddSeriesAttribute(self, groupElem):
        """
        Add the specified attribute to the set of series attributes to be scanned.

        Parameters
        ----------
        groupElem : (int,int) tuple
            (group,elem) tag for the attribute to be added.
        """
        self.seriesTags.append(groupElem)

    def AddAcquisitionAttribute(self, groupElem):
        """
        Add the specified attribute to the set of acquisition attributes to be scanned.

        Parameters
        ----------
        groupElem : (int,int) tuple
            (group,elem) tag for the attribute to be added.
        """
        self.acquisitionTags.append(groupElem)

    class ScanHelper(object):
        """
        Callable passed to the DicomScan function.
        """

        def __init__(self, patientTags, studyTags, seriesTags, acquisitionTags):
            """
            Initialize the helper.
            """

            # Copy references to tag lists
            self.patientTags     = patientTags
            self.studyTags       = studyTags
            self.seriesTags      = seriesTags
            self.acquisitionTags = acquisitionTags

            # Initialize the results
            self.results = {}

        def __call__(self, path, tvDict):
            """
            Function called once for each file to update results.
            """

            # Skip if the unique IDs are not present
            if (StudyInstanceUIDTag  not in tvDict or
                SeriesInstanceUIDTag  not in tvDict or
                SOPInstanceUIDTag not in tvDict): return

            #-------------------------------------------------

            # Get the patient ID, or None if it is absent
            patId = tvDict[PatientIDTag] if PatientIDTag in tvDict else None

            # If this is the first time we've seen this patient
            if patId not in self.results:

                # Create a new entry for it
                self.results[patId] = [{},{}]
                pat = self.results[patId]

                # Fill it with the requested attributes
                for ge in self.patientTags: pat[0][ge] = tvDict[ge] if ge in tvDict else None

            else: pat = self.results[patId]

            #-------------------------------------------------

            # Get the study UID
            stdUid = tvDict[StudyInstanceUIDTag]

            # If this is the first time we've seen this study
            if stdUid not in pat[1]:

                # Create a new entry for it
                pat[1][stdUid] = [{},{}]
                std = pat[1][stdUid]

                # Fill it with the requested attributes
                for ge in self.studyTags: std[0][ge] = tvDict[ge] if ge in tvDict else None

            else: std = pat[1][stdUid]

            #-------------------------------------------------

            # Get the series UID
            serUid = tvDict[SeriesInstanceUIDTag]

            # If this is the first time we've seen this series
            if serUid not in std[1]:

                # Create a new entry for it
                std[1][serUid] = [{},{}]
                ser = std[1][serUid]

                # Fill it with the requested attributes
                for ge in self.seriesTags: ser[0][ge] = tvDict[ge] if ge in tvDict else None

            else: ser = std[1][serUid]

            #-------------------------------------------------

            # Get the instance UID
            instUid = tvDict[SOPInstanceUIDTag]

            # Skip if we've already processed this instance (duplicate)
            if instUid in ser[1]: return

            # Create the instance
            ser[1][instUid] = [{},path]
            inst = ser[1][instUid]

            # Set the instance level tags
            instTags = self.acquisitionTags.copy()
            instTags.append(SOPClassUIDTag)
            instTags.append(ImageRowsTag)
            instTags.append(ImageColumnsTag)
            instTags.append(PixelSpacingTag)
            instTags.append(InstanceNumberTag)
            instTags.append(SlicePositionTag)
            instTags.append(ImagePositionPatientTag)
            instTags.append(ImageOrientationPatientTag)
            instTags.append(EchoTimeTag)

            # Fill it with the requested attributes
            for ge in instTags: inst[0][ge] = tvDict[ge] if ge in tvDict else None

    def Scan(self, paths):
        """
        Execute the scan on the specified paths.

        Parameters
        ----------
        paths : str or sequence of str
            Names of files and directories to be scanned. Can be a
            mixture of both.

        Returns
        -------
        The raw results of the scan in a patient/study/series/instance hierarchy.
        This is typically discarded, unless many per-instance attributes are needed.

        Note that the returned results are different from the patient/study/series/
        acquisition hierarchy that is held in the object following the scan, and
        saved by the Save method. That data saved by the Save method retains very
        little per-instance data, only the Image Position (Patient), SOP Instance
        UID, and path for each instance. In contrast, the returned results hold all
        data specified to be gathered at the acquisition level for every instance.

        The returned object is a dict from Patient ID to a list containing two dicts.
        The first maps (group,element) pairs to their values for attributes collected
        at the patient level.

        The second maps Study Instance UIDs to a list containing two dicts. The first
        maps (group,element) pairs to their values for attributes collected at the
        study level.

        The second maps Series Instance UIDs to a list containing two dicts. The first
        maps (group,element) pairs to their values for attributes collected at the
        series level.

        The second maps SOP Instance UIDs to a list containing a dict and a string. The
        dict maps (group,element) pairs to their values for attributes collected at the
        acquisition level for each instance. The string is the path to the instance.
        """

        # Make sure the input is a path or a list of paths
        if not isinstance(paths, str) and (not isinstance(paths, collections.
            Sequence) or len([p for p in paths if not isinstance(p, str)]) != 0):
            raise TypeError("Expected a path or a sequence of paths")

        # Is the argument a group/element tuple?
        def isGroupElem(ge): return (isinstance(ge, tuple) and len(ge) == 2
            and isinstance(ge[0], int) and isinstance(ge[1], int))

        # Make sure the inputs are lists of group/elem tuples
        if (not isinstance(self.patientTags, collections.Sequence) or len(
            [ge for ge in self.patientTags if not isGroupElem(ge)]) != 0):
            raise TypeError("Expected a list of group/elem tuples in patient tags")

        if (not isinstance(self.studyTags, collections.Sequence) or len(
            [ge for ge in self.studyTags if not isGroupElem(ge)]) != 0):
            raise TypeError("Expected a list of group/elem tuples in study tags")

        if (not isinstance(self.seriesTags, collections.Sequence) or len(
            [ge for ge in self.seriesTags if not isGroupElem(ge)]) != 0):
            raise TypeError("Expected a list of group/elem tuples in series tags")

        if (not isinstance(self.acquisitionTags, collections.Sequence) or len(
            [ge for ge in self.acquisitionTags if not isGroupElem(ge)]) != 0):
            raise TypeError("Expected a list of group/elem tuples in acquisition tags")

        # The maximum amount by which unit vectors can deviate
        # from having unit norm
        unitVecMagDev = 0.02

        # The distance netween adjacent slices must not be less
        # than this amount
        minAllowedSliceGap = 0.02; # mm

        # For any three adjacent slices, the difference between the
        # gaps separating the first two and the second two must not
        # be greater than this amount
        maxAllowedGapDev = 0.02; # mm

        # The direction cosine between vectors separating adjacent
        # pairs of slices and between these vectors and the acquisition
        # axis must not differ from 1 by more than this amount
        maxAllowedDirCosDev = 0.0075; # 7.0 degrees

        # Assemble the full list of attributes to be queried
        groupElem = []
        groupElem.extend(self.patientTags)
        groupElem.extend(self.studyTags)
        groupElem.extend(self.seriesTags)
        groupElem.extend(self.acquisitionTags)
        groupElem.extend([
            PatientIDTag,
            StudyInstanceUIDTag,
            SeriesInstanceUIDTag,
            SOPClassUIDTag,
            SOPInstanceUIDTag,
            ImageRowsTag,
            ImageColumnsTag,
            PixelSpacingTag,
            InstanceNumberTag,
            SlicePositionTag,
            ImagePositionPatientTag,
            ImageOrientationPatientTag,
            EchoTimeTag])

        # Create the helper callable
        scanHelp = DicomScanner.ScanHelper(self.patientTags, self.studyTags, self.seriesTags, self.acquisitionTags)

        # Call the scanner
        DicomScan(paths, groupElem, scanHelp)

        # Save the results
        results = scanHelp.results

        # Clear the patients
        self.patients = []

        # Get the sorted indices of patient IDs
        patList = [(k,v) for (k,v) in results.items()]
        patKeys = [kv[0] if kv[0] else "None" for kv in patList]
        patIdxs = sorted(range(len(patKeys)), key=patKeys.__getitem__)

        # Assemble the sorted list of patients
        for patIdx in patIdxs:

            # Get a shortcut to the patient
            pat = patList[patIdx]

            # Initialize the new patient with the patient ID
            newPat = DicomPatient()
            newPat.attrs[PatientIDTag] = pat[0]

            # Add in all the other patient attributes
            for ge in self.patientTags: newPat.attrs[ge] = pat[1][0][ge]

            # Get the sorted indices of study dates
            stdDtTag = StudyDateTag
            stdList =  [(k,v) for (k,v) in pat[1][1].items()]
            stdKeys = [kv[1][0][stdDtTag] if kv[1][0][stdDtTag] else "None" for kv in stdList]
            stdIdxs = sorted(range(len(stdKeys)), key=stdKeys.__getitem__, reverse=True)

            # Assemble the sorted list of studies
            for stdIdx in stdIdxs:

                # Get a shortcut to the study
                std = stdList[stdIdx]

                # Initialize the new study with the study UID
                newStd = DicomStudy()
                newStd.attrs[StudyInstanceUIDTag] = std[0]

                # Add in all the other study attributes
                for ge in self.studyTags: newStd.attrs[ge] = std[1][0][ge]

                # Get the series
                serList =  [(k,v) for (k,v) in std[1][1].items()]

                # Initialize the list of series data. For each series, we want to collect
                # all the different acquisitions. For each acquisition, we want the series
                # UID, a representative instance UID, the dimensions, the spacing, whether
                # the spacing is homogeneous, the acquisition axis, and the list of file paths
                serData = []

                # First pass over series to get acquisitions
                for ser in serList:

                    # Get the modality
                    mod = ser[1][0][SeriesModalityTag]

                    # Acquisitions for this series
                    acqData = []

                    # Collect the instance attributes
                    rows = []; cols = []; pix = []; typ  = []
                    prot = []; inn  = []; slp = []; ipp  = []
                    iop  = []; adt  = []; atm = []; inst = []; paths = []
                    et   = []
                    for kv in ser[1][1].items():
                        rows .append(kv[1][0][ImageRowsTag])
                        cols .append(kv[1][0][ImageColumnsTag])
                        pix  .append(kv[1][0][PixelSpacingTag])
                        typ  .append(kv[1][0][ImageTypeTag])
                        prot .append(kv[1][0][ProtocolNameTag])
                        inn  .append(kv[1][0][InstanceNumberTag])
                        slp  .append(kv[1][0][SlicePositionTag])
                        ipp  .append(kv[1][0][ImagePositionPatientTag])
                        iop  .append(kv[1][0][ImageOrientationPatientTag])
                        adt  .append(kv[1][0][AcquisitionDateTag])
                        atm  .append(kv[1][0][AcquisitionTimeTag])
                        et   .append(kv[1][0][EchoTimeTag])
                        inst .append(kv[0])
                        paths.append(kv[1][1])

                    # Get the number of files
                    nFiles = len(paths)

                    # Return a single acquisition vector
                    def GetAcquisitionVector(singleIop):
                        if singleIop:
                            singleIopNp = np.array(singleIop)
                            singleAcqVec = list(np.cross(singleIopNp[:3], singleIopNp[3:]))
                            singleAcqVec = singleAcqVec[::-1] # Convert to Z,Y,X
                        else: singleAcqVec = None
                        return singleAcqVec

                    # If there is more than one file
                    if nFiles > 1:

                        # Remove small rounding errors for pixel spacing
                        # and image orientation (patient) for hashing
                        pix = [[round(p1,4) for p1 in p2] if p2 else p2 for p2 in pix]
                        iop = [[round(p1,4) for p1 in p2] if p2 else p2 for p2 in iop]
                        et  = [round(float(e),4) if e is not None and e != '' else None for e in et]

                        # Collect the attributes together for allocation to acquisitions
                        instAttrs = list(zip(rows, cols, pix, typ, prot, inn, slp, ipp, iop, adt, atm, inst, paths, et))

                        # Collect instance attributes into unique acquisition groups
                        uniqAttrs = {}
                        for i in instAttrs:

                            # Construct the hash key
                            hsh = (i[0], i[1],
                                tuple(i[2]) if i[2] else None,
                                tuple(i[3]) if i[3] else None, i[4],
                                tuple(i[8]) if i[8] else None, i[13])

                            # Create a new entry if we don't have it
                            if hsh in uniqAttrs: uniqAttrs[hsh].append(i)
                            else: uniqAttrs[hsh] = [i]

                        # Loop over the acquisitions
                        for u in uniqAttrs.values():

                            # Get instance attributes to the current unique set
                            (rows2, cols2, pix2, typ2, prot2, inn2, slp2, ipp2,
                                iop2, adt2, atm2, inst2, paths2, et2) = list(zip(*u))

                            # Does this sub-series have more than one instance?
                            if len(inst2) > 1:

                                # Figure out what kind of spacing data we have
                                allImagePos = (all(ipp2) and all(iop2))
                                allAcqDAndT = (all(adt2) and all(atm2))
                                allSlicePos = all(slp2)
                                allInstNums = all(inn2)

                                # If we have image orientations and positions,
                                # check if the orientations are also consistent
                                if allImagePos:

                                    # Get the patient position/orientation data
                                    ipp2Np = np.array(ipp2)
                                    iop2Np = np.array(iop2)

                                    # Compute the acquisition vectors
                                    acqVecs = np.cross(iop2Np[:,:3], iop2Np[:,3:])

                                    # Update image orientations and positions status based on homogeneity
                                    allImagePos = (all(np.fabs(np.linalg.norm(iop2Np[:,:3], axis=1) - 1.0) <= unitVecMagDev) and
                                                   all(np.fabs(np.linalg.norm(iop2Np[:,3:], axis=1) - 1.0) <= unitVecMagDev) and
                                                   all(np.fabs(np.linalg.norm(acqVecs, axis=1) - 1.0) <= unitVecMagDev) and
                                                   all(np.fabs(np.dot(acqVecs[0,:], acqVecs.T) - 1.0) <= maxAllowedDirCosDev))

                                # If we have 3D positions and orientations
                                if allImagePos:

                                    # Compute the positions along the acquisition axis
                                    ipp2Np -= ipp2Np[0,:] * np.ones((ipp2Np.shape[0],1), dtype=np.float32)
                                    slicePos = np.dot(ipp2Np, acqVecs[0,:].T)

                                elif allSlicePos: # Just scalar positions

                                    # Get the slice positions
                                    slicePos = np.array(slp2, np.float32)

                                elif allInstNums: # Just image numbers

                                    # Get the image numbers
                                    slicePos = np.array(inn2, np.float32)

                                # Do we have more than two instances?
                                mtt = len(inst2) > 2

                                # Initialize defaults
                                nPosTrpls = nDTPTrpls = nINNTrpls = -1
                                posIdxs   = dtpIdxs   = innIdxs   = None
                                posGaps   = dtpGaps   = innGaps   = None
                                posTrpls  = dtpTrpls  = innTrpls  = None

                                # Get the ordering data for sorting via positions w/wo acq date/
                                # times vs sorting by instance numbers. We will consider all
                                if allImagePos or allSlicePos:
                                    posIdxs = np.argsort(slicePos)
                                    posGaps = np.diff(slicePos[posIdxs])
                                    if mtt:
                                        posTrpls = np.abs(np.diff(posGaps))
                                        posTrpls = np.less(posTrpls, maxAllowedGapDev)
                                        posTrpls = np.logical_and(posTrpls, np.greater(
                                            np.abs(posGaps[:-1]), minAllowedSliceGap))
                                        nPosTrpls = np.sum(posTrpls)
                                if allImagePos and allAcqDAndT:
                                    dtpKeys = [(ad,at,sp) for ad, at, sp in zip(adt2,atm2,slicePos)]
                                    dtpIdxs = np.array(sorted(range(len(dtpKeys)), key=dtpKeys.__getitem__))
                                    dtpGaps = np.diff(slicePos[dtpIdxs])
                                    if mtt:
                                        dtpTrpls = np.abs(np.diff(dtpGaps))
                                        dtpTrpls = np.less(dtpTrpls, maxAllowedGapDev)
                                        dtpTrpls = np.logical_and(dtpTrpls, np.greater(
                                            np.abs(dtpGaps[:-1]), minAllowedSliceGap))
                                        nDTPTrpls = np.sum(dtpTrpls)
                                if allInstNums:
                                    innIdxs = np.argsort(inn2)
                                    innGaps = np.diff(slicePos[innIdxs])
                                    if mtt:
                                        innTrpls = np.abs(np.diff(innGaps))
                                        innTrpls = np.less(innTrpls, maxAllowedGapDev)
                                        innTrpls = np.logical_and(innTrpls, np.greater(
                                            np.fabs(innGaps[:-1]), minAllowedSliceGap))
                                        nINNTrpls = np.sum(innTrpls)

                                # Initially assume no spacing info
                                gaps = trpls = None

                                # Figure out the ordering to use and get the triples
                                if mtt and (allImagePos or allSlicePos or allInstNums):

                                    # Choose the ordering information that
                                    # provides the least failed triples
                                    sumTrpls = [nPosTrpls,nDTPTrpls,nINNTrpls]
                                    bestIdx  = sorted(range(len(sumTrpls)), key=sumTrpls.__getitem__)[-1]

                                    # Set best idxs, gaps, triples
                                    allIdxs  = [posIdxs,dtpIdxs,innIdxs]
                                    allGaps  = [posGaps,dtpGaps,innGaps]
                                    allTrpls = [posTrpls,dtpTrpls,innTrpls]
                                    idxs = allIdxs[bestIdx]
                                    if allTrpls[bestIdx] is not None and np.sum(
                                        allTrpls[bestIdx]) > len(inst2) / 2:
                                        gaps = allGaps[bestIdx]
                                        trpls = allTrpls[bestIdx]

                                else:

                                    # Create a dummy ordering
                                    idxs = [i for i in range(len(inst2))]

                                # Re-order the positions, re-ordering to Z,Y,X order
                                if allImagePos: pos2 = [ipp2[i][::-1] for i in idxs]
                                elif allSlicePos: pos2 = [slp2[i] for i in idxs]
                                else: pos2 = [ipp2[i] for i in idxs]

                                # Re-order the instance UIDs, paths
                                inst2 = [inst2[i] for i in idxs]
                                paths2 = [paths2[i] for i in idxs]

                                # If there is no inhomogeneity
                                if trpls is None or np.all(trpls):

                                    # Set the spacing data
                                    if gaps is not None:
                                        spcInfo = np.median(gaps)
                                        if pix2[0]: spcInfo = [spcInfo, *pix2[0]]
                                        else: spcInfo = [spcInfo, None, None]
                                    else:
                                        if pix2[0]: spcInfo = [None, *pix2[0]]
                                        else: spcInfo = None

                                    # Reverse direction if necessary
                                    if (spcInfo is not None and
                                        spcInfo[0] is not None and
                                        spcInfo[0] < 0.0):
                                        spcInfo[0] = -spcInfo[0]
                                        pos2 = pos2[::-1]
                                        inst2 = inst2[::-1]
                                        paths2 = paths2[::-1]

                                    # Record the single acquisition
                                    acqData.append([])
                                    acqData[-1].append(ser[0])
                                    acqData[-1].append(inst2[0])
                                    acqData[-1].append([len(inst2),rows2[0],cols2[0]])
                                    acqData[-1].append(spcInfo)
                                    acqData[-1].append(gaps is not None)
                                    acqData[-1].append(GetAcquisitionVector(iop2[0]))
                                    acqData[-1].append(pos2)
                                    acqData[-1].append(inst2)
                                    acqData[-1].append(paths2)
                                    acqData[-1].append(et2)

                                else: # Need to split

                                    # We initially start outside an acquisition
                                    inAcq = False

                                    # Initially set all slices unused
                                    unused = np.full((len(inst2),), True, np.bool)

                                    # Loop over the triples
                                    for t in range(len(trpls)):

                                        # Are we not in an acquisition?
                                        if not inAcq:

                                            # Is the current triple valid?
                                            if trpls[t]:

                                                # Starting a new acquisition
                                                pos3 = [pos2[t],pos2[t+1],pos2[t+2]]
                                                inst3 = [inst2[t],inst2[t+1],inst2[t+2]]
                                                paths3 = [paths2[t],paths2[t+1],paths2[t+2]]
                                                unused[t] = unused[t+1] = unused[t+2] = False
                                                inAcq = True

                                        else: # In an acquisition

                                            # Is the current triple valid?
                                            if trpls[t]:

                                                # Add current slice to acquisition
                                                pos3.append(pos2[t+2])
                                                inst3.append(inst2[t+2])
                                                paths3.append(paths2[t+2])
                                                unused[t+2] = False

                                            # If invalid or last triple
                                            if not trpls[t] or t == len(trpls)-1:

                                                # Set the spacing data
                                                spcInfo = gaps[t]
                                                if pix2[0]: spcInfo = [spcInfo, *pix2[0]]
                                                else: spcInfo = [spcInfo, None, None]

                                                # Reverse direction if necessary
                                                if spcInfo[0] < 0.0:
                                                    spcInfo[0] = -spcInfo[0]
                                                    pos3 = pos3[::-1]
                                                    inst3 = inst3[::-1]
                                                    paths3 = paths3[::-1]

                                                # Record the single acquisition
                                                acqData.append([])
                                                acqData[-1].append(ser[0])
                                                acqData[-1].append(inst3[0])
                                                acqData[-1].append([len(inst3),rows2[0],cols2[0]])
                                                acqData[-1].append(spcInfo)
                                                acqData[-1].append(True)
                                                acqData[-1].append(GetAcquisitionVector(iop2[0]))
                                                acqData[-1].append(pos3)
                                                acqData[-1].append(inst3)
                                                acqData[-1].append(paths3)

                                                # End of acquisition
                                                inAcq = False

                                    # Check for unused slices
                                    if np.any(unused):

                                        # Collect them into a runt acquisition
                                        acqData.append([])
                                        acqData[-1].append(ser[0])
                                        acqData[-1].append(inst2[0])
                                        acqData[-1].append([np.sum(unused),rows2[0],cols2[0]])
                                        acqData[-1].append([None, *pix2[0]] if pix2[0] else None)
                                        acqData[-1].append(False)
                                        acqData[-1].append(GetAcquisitionVector(iop2[0]))
                                        acqData[-1].append([i[0] for i in zip(ipp2, unused) if i[1]])
                                        acqData[-1].append([i[0] for i in zip(inst2, unused) if i[1]])
                                        acqData[-1].append([i[0] for i in zip(paths2, unused) if i[1]])

                            else:

                                # Just a single instance
                                acqData.append([])
                                acqData[-1].append(ser[0])
                                acqData[-1].append(inst2[0])
                                acqData[-1].append([rows2[0], cols2[0]] if rows2[0] and cols2[0] else None)
                                acqData[-1].append(list(pix2[0]) if pix2[0] else None)
                                acqData[-1].append(None)
                                acqData[-1].append(GetAcquisitionVector(iop2[0]))
                                acqData[-1].append(ipp2[::-1] if ipp2 else slp2)
                                acqData[-1].append(inst2)
                                acqData[-1].append(paths2)

                    else: # Just a single file

                        # TODO: Add branches here for mammo, tomo, spect

                        # Only one file
                        acqData.append([])
                        acqData[-1].append(ser[0])
                        acqData[-1].append(inst[0])
                        acqData[-1].append([rows[0], cols[0]] if rows[0] and cols[0] else None)
                        acqData[-1].append(list(pix[0]) if pix[0] else None)
                        acqData[-1].append(None)
                        acqData[-1].append(GetAcquisitionVector(iop[0]))
                        acqData[-1].append(ipp[::-1] if ipp else slp)
                        acqData[-1].append(inst)
                        acqData[-1].append(paths)
                        acqData[-1].append(et)

                    # Get the sorted indices of acquisitions
                    acqKeys = []
                    for a in acqData:
                        dc = ser[1][1][a[1]][0]
                        da = dc[AcquisitionDateTag] if AcquisitionDateTag in dc else None
                        tm = dc[AcquisitionTimeTag] if AcquisitionTimeTag in dc else None
                        acqKeys.append((da if da else "19000101", tm if tm else "000000"))
                    acqIdxs = sorted(range(len(acqKeys)), key=acqKeys.__getitem__)
                    acqData = [acqData[ai] for ai in acqIdxs]

                    # Save the results on this series
                    serData.append((mod,acqData))

                # Get the sorted indices of series
                serKeys = []
                for s in serData:
                    dc = std[1][1][s[1][0][0]][1][s[1][0][1]][0]
                    da = dc[AcquisitionDateTag] if AcquisitionDateTag in dc else None
                    tm = dc[AcquisitionTimeTag] if AcquisitionTimeTag in dc else None
                    serKeys.append((s[0] ,da if da else "19000101", tm if tm else "000000"))
                    serIdxs = sorted(range(len(serKeys)), key=serKeys.__getitem__)

                # Assemble the sorted list of series
                for serIdx in serIdxs:

                    # Get a shortcut to the series
                    serL = serList[serIdx]
                    serD = serData[serIdx][1]

                    # Initialize the new series with the series UID
                    newSer = DicomSeries()
                    newSer.attrs[SeriesInstanceUIDTag] = serL[0]

                    # Add in all the other series attributes
                    for ge in self.seriesTags: newSer.attrs[ge] = serL[1][0][ge]

                    # Assemble the sorted list of acquisitions
                    for acq in serD:

                        # Initialize the new acquisition
                        newAcq = DicomAcquisition()

                        # Get the representative instance
                        inst = serL[1][1][acq[1]]

                        # Add in all the acquisition attributes
                        for ge in self.acquisitionTags: newAcq.attrs[ge] = inst[0][ge]

                        # Append the additional acquisition data
                        newAcq.dimensions = acq[2]
                        newAcq.spacing = acq[3]
                        newAcq.homogeneous = acq[4]
                        newAcq.acquisitionAxis = acq[5]
                        newAcq.imgPosPats = acq[6]
                        newAcq.sopInstUIDs = acq[7]
                        newAcq.children = acq[8]

                        # Set the parent pointer and store the assembled acquisition
                        newAcq.parent = wr.ref(newSer)
                        newSer.children.append(newAcq)

                    # Set the parent pointer and store the assembled series
                    newSer.parent = wr.ref(newStd)
                    newStd.children.append(newSer)

                # Set the parent pointer and store the assembled study
                newStd.parent = wr.ref(newPat)
                newPat.children.append(newStd)

            # Store the assembled patient
            self.patients.append(newPat)

        # Return the raw results
        return results

    def Rescan(self):
        """
        Rescan all the files to update internal data structures.
        """

        # Initialize the files to scan
        filesToScan = []

        # Collect all the files
        for pat in self.patients:
            for std in pat:
                for ser in std:
                    for acq in ser:
                        filesToScan.extend(acq.GetFilePaths())

        # Execute the scan
        self.Scan(filesToScan)

    def ConvertPaths(self, oldPart, newPart, oldSep="", newSep=""):
        """
        Replace all instances of oldPart with newPart
        and oldSep with newSep in all file paths.
        """

        # Loop over all acquisitions
        for pat in self.patients:
            for std in pat:
                for ser in std:
                    for acq in ser:

                        # Replace old paths with converted paths
                        acq.children = [fp.replace(oldPart, newPart).replace(
                            oldSep, newSep) for fp in acq.children]

    def GetPatients(self):
        """
        Return the list of DicomPatient objects at the top of the
        patient/study/series/acquisition hierarchy. Alternately,
        iterate over self if you don't need the whole list.
        """

        # Return the list of patients
        return self.patients

    def Print(self, file=sys.stdout):
        """
        Print the the patient/study/series/acquisition hierarchy to file.

        Parameters
        ----------
        file : str or file handle (optional)
            Name of the file, or handle to the file, to which the results should
            be written. If not specified, results are printed to sys.stdout.
        """

        # Convert file name to handle if necessary
        file = open(file, "w") if isinstance(file, str) else file

        # Read dictionary to get attribute names
        globInst = gdcm.Global.GetInstance()
        globInst.LoadResourcesFiles()
        dicts = globInst.GetDicts()

        # Set the indent
        ind = "    "

        # Loop over patients
        patNo = 0
        for pat in self.patients:

            # Print out requested attributes
            patNo += 1
            patInd = ind
            print("PATIENT {} of {}".format(patNo, len(self.patients)), file=file)
            for ge in [PatientIDTag,*self.patientTags]:
                if(pat.attrs[ge]): print("{}{}: {}".format(patInd, dicts.GetDictEntry(
                    gdcm.Tag(ge[0], ge[1])).GetName(), pat.attrs[ge]), file=file)

            # Loop over studies
            stdNo = 0
            for std in pat:

                # Print out requested attributes
                stdNo += 1
                stdInd = "{}{}".format(patInd, ind)
                print("{}STUDY {} of {}".format(patInd, stdNo, len(pat.GetChildren())), file=file)
                for ge in [StudyInstanceUIDTag,*self.studyTags]:
                    if(std.attrs[ge]): print("{}{}: {}".format(stdInd, dicts.GetDictEntry(
                        gdcm.Tag(ge[0], ge[1])).GetName(), std.attrs[ge]), file=file)

                # Loop over series
                serNo = 0
                for ser in std:

                    # Print out requested attributes
                    serNo += 1
                    serInd = "{}{}".format(stdInd, ind)
                    print("{}SERIES {} of {}".format(stdInd, serNo, len(std.GetChildren())), file=file)
                    for ge in [SeriesInstanceUIDTag,*self.seriesTags]:
                        if(ser.attrs[ge]): print("{}{}: {}".format(serInd, dicts.GetDictEntry(
                            gdcm.Tag(ge[0], ge[1])).GetName(), ser.attrs[ge]), file=file)

                    # Loop over acquisitions
                    acqNo = 0
                    for acq in ser:

                        # Skip if nothing to print out
                        acqTags = [acq.attrs[ge] for ge in self.acquisitionTags]
                        if not any(acqTags) and \
                           not acq.dimensions and \
                           not acq.spacing and \
                           not acq.homogeneous and \
                           not acq.acquisitionAxis: continue

                        # Print out requested attributes
                        acqNo += 1
                        acqInd = "{}{}".format(serInd, ind)
                        print("{}ACQUISITION {} of {}".format(serInd, acqNo, len(ser.GetChildren())), file=file)
                        for ge in self.acquisitionTags:
                            if(acq.attrs[ge]): print("{}{}: {}".format(acqInd, dicts.GetDictEntry(
                                gdcm.Tag(ge[0], ge[1])).GetName(), acq.attrs[ge]), file=file)

                        # Print out dimensions, spacing, homogeneity, and acqusition axis
                        if acq.GetDimensions(): print("{}Dimensions: {}".format(acqInd, acq.GetDimensions(Formats.String)), file=file)
                        if acq.GetSpacing(): print("{}Spacing: {}".format(acqInd, acq.GetSpacing(Formats.String)), file=file)
                        if acq.GetHomogeneous() is not None: print("{}Homogeneous: {}".format(acqInd, acq.GetHomogeneous(Formats.String)), file=file)
                        if acq.GetAcquisitionAxis(): print("{}Acquisition Axis: {}".format(acqInd, acq.GetAcquisitionAxis(Formats.String)), file=file)

    def Save(self, path):
        """
        Save the scan results to a file on disk.

        Parameters
        ----------
        path : str
            Path to which output should be written.
        """

        # Make sure the input is a path
        if not isinstance(path, str): raise TypeError(
            "Expected a string containing a path")

        # Remove parent weak references
        for pat in self.patients:
            for std in pat:
                std.parent = None
                for ser in std:
                    ser.parent = None
                    for acq in ser:
                        acq.parent = None

        # Set the save version
        version = 1

        # Aggregate the instance attributes to save
        outData = [version,self.patientTags,self.studyTags,
            self.seriesTags,self.acquisitionTags,self.patients]

        # Save the results
        with open(path, "wb") as output:
            pickle.dump(outData, output, -1)

        # Restore parent weak references,
        # as pickling seems to destroy them
        for pat in self.patients:
            for std in pat:
                std.parent = wr.ref(pat)
                for ser in std:
                    ser.parent = wr.ref(std)
                    for acq in ser:
                        acq.parent = wr.ref(ser)

    def Load(self, path):
        """
        Load the scan results from a file on disk.

        Parameters
        ----------
        path : str
            Path from which output should be read.
        """

        # Make sure the input is a path
        if not isinstance(path, str): raise TypeError(
            "Expected a string containing a path")

        # Read the results
        with open(path, "rb") as input:
            inData = pickle.load(input)

        # Make sure we can read it
        if len(inData) != 6: raise TypeError(
            "Unexpected data format")

        # Restore the instance attributes
        version              = inData[0]
        self.patientTags     = inData[1]
        self.studyTags       = inData[2]
        self.seriesTags      = inData[3]
        self.acquisitionTags = inData[4]
        self.patients        = inData[5]

        # Restore parent weak references
        for pat in self.patients:
            for std in pat:
                std.parent = wr.ref(pat)
                for ser in std:
                    ser.parent = wr.ref(std)
                    for acq in ser:
                        acq.parent = wr.ref(ser)
