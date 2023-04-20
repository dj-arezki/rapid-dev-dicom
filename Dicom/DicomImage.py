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
Class for reading and querying DICOM images.

"""

import numpy as np
from numpy.compat import asbytes

import gdcm

from Dicom.DicomFile import DicomFile

class DicomImage(DicomFile):
    """
    Class for reading and querying DICOM images.    
    """
    
    def __init__(self, fileName=None):
        """
        Initialize the object, and optionally read an image from file.
 
        Parameters
        ----------
        fileName : str
          The (optional) name of the file to be read.
        """
        
        # Call the base class function
        super(DicomImage, self).__init__(fileName)

        # Initialize image attribute
        self.image = None
         
    def GetReader(self):
        """
        Polymorphic function for returning a kind of reader.
        """
        return gdcm.ImageReader()

    def GetPixelData(self, native = False):
        """
        Get a numpy array with the image pixel data.
 
        Parameters
        ----------
        native : bool
          Should the pixel data be returned in its native format? If False
              (the default), the pixel data will be returned as numpy.float32,
              after applying the rescale slope and intercept if they are
              present.
        """
        
        # Make sure we've read something
        if self.reader is None:
            raise RuntimeError('No DICOM file has been loaded')
            
        # Get the image if not done yet
        if self.image is None:
            self.image = self.reader.GetImage()
            
        # Figure out the type of the data
        gdcmTypemap = {
            gdcm.PixelFormat.INT8:     np.int8,
            gdcm.PixelFormat.UINT8:    np.uint8,
            gdcm.PixelFormat.UINT16:   np.uint16,
            gdcm.PixelFormat.INT16:    np.int16,
            gdcm.PixelFormat.UINT32:   np.uint32,
            gdcm.PixelFormat.INT32:    np.int32,
            gdcm.PixelFormat.FLOAT32:  np.float32,
            gdcm.PixelFormat.FLOAT64:  np.float64
        }
        pixelFormat = self.image.GetPixelFormat().GetScalarType()
        if pixelFormat not in gdcmTypemap:
            raise RuntimeError('DICOM image format not supported')
            
        # Get the raw data
        data = np.frombuffer(self.image.GetBuffer().encode("utf-8",
            errors="surrogateescape"), dtype=gdcmTypemap[pixelFormat])
        
        # Reshape it
        self.shape = self.image.GetDimensions()
        data = data.reshape(self.shape[::-1])
        
        # Return to caller as appropriate type
        if not native: 
            
            # Try to get the rescaling
            rescaleIcept = self.GetAttribute((0x0028,0x1052))
            rescaleSlope = self.GetAttribute((0x0028,0x1053))

            rescaleIcept = float(rescaleIcept) if rescaleIcept else 0.0
            rescaleSlope = float(rescaleSlope) if rescaleSlope else 1.0

            return rescaleSlope * data.astype(np.float32) + rescaleIcept

        else: return data
