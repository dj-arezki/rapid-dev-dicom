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
Class for reading and querying DICOM data sets.

"""

import gdcm

from Dicom.DicomData import DicomData

class DicomFile(DicomData):
    """
    Class for reading and querying DICOM data sets.    
    """
    
    def __init__(self, fileName=None):
        """
        Initialize the object, and optionally read a data set from file.
 
        Parameters
        ----------
        fileName : str
          The (optional) name of the file to be read.
        """

        # Call the base class function
        super(DicomFile, self).__init__()
        
        # Read file if file name specified
        if fileName is not None: self.Read(fileName)
         
    def GetReader(self):
        """
        Polymorphic function for returning a kind of reader.
        """
        return gdcm.Reader()

    def Read(self, fileName):
        """
        Read a DICOM file.
 
        Parameters
        ----------
        fileName : str
          The name of the file to be read.
        """

        # Get the polymorphic reader
        self.reader = self.GetReader()
        
        # Try to read the DICOM file
        self.reader.SetFileName(fileName)
        if not self.reader.Read():
            raise IOError('Failed to read DICOM file {}'.format(fileName))
            
        # Extract all the data for later use
        self.file = self.reader.GetFile()
        self.dataSet = self.file.GetDataSet()
         
        # For parsing attributes
        self.strFilter = gdcm.StringFilter()
        self.strFilter.SetFile(self.file)
