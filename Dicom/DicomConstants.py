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
Attributes for the DICOM hierarchy.
"""

#------------------
# Common attributes
#------------------

TransferSyntaxUIDTag       = (0x0002,0x0010)
ImageTypeTag               = (0x0008,0x0008)
SOPClassUIDTag             = (0x0008,0x0016)
SOPInstanceUIDTag          = (0x0008,0x0018)
StudyDateTag               = (0x0008,0x0020)
AcquisitionDateTag         = (0x0008,0x0022)
AcquisitionTimeTag         = (0x0008,0x0032)
AccessionNumberTag         = (0x0008,0x0050)
SeriesModalityTag          = (0x0008,0x0060)
ManufacturerTag            = (0x0008,0x0070)
InstitutionNameTag         = (0x0008,0x0080)
StudyDescriptionTag        = (0x0008,0x1030)
SeriesDescriptionTag       = (0x0008,0x103E)
ManufacturerModelTag       = (0x0008,0x1090)
PatientsNameTag            = (0x0010,0x0010)
PatientIDTag               = (0x0010,0x0020)
PatientsSexTag             = (0x0010,0x0040)
PatientsAgeTag             = (0x0010,0x1010)
PatientCommentsTag         = (0x0010,0x4000)
SliceThicknessTag          = (0x0018,0x0050)
KVPTag                     = (0x0018,0x0060)
EchoTimeTag                = (0x0018,0x0081)
SoftwareVersionsTag        = (0x0018,0x1020)
ProtocolNameTag            = (0x0018,0x1030)
ContrastBolusAgentTag      = (0x0018,0x1048)
GantryTiltTag              = (0x0018,0x1120)
XRayTubeCurrentTag         = (0x0018,0x1151)
ExposureTag                = (0x0018,0x1152)
FilterTypeTag              = (0x0018,0x1160)
ConvolutionKernelTag       = (0x0018,0x1210)
MultiEnergyCTAcqTag        = (0x0018,0x9361)
StudyInstanceUIDTag        = (0x0020,0x000D)
SeriesInstanceUIDTag       = (0x0020,0x000E)
InstanceNumberTag          = (0x0020,0x0013)
ImagePositionPatientTag    = (0x0020,0x0032)
ImageOrientationPatientTag = (0x0020,0x0037)
FrameOfReferenceUIDTag     = (0x0020,0x0052)
SlicePositionTag           = (0x0020,0x1041)
ImageCommentsTag           = (0x0020,0x4000)
ImageRowsTag               = (0x0028,0x0010)
ImageColumnsTag            = (0x0028,0x0011)
PixelSpacingTag            = (0x0028,0x0030)
BurnedInAnnotationsTag     = (0x0028,0x0301)
LossyImageCompressionTag   = (0x0028,0x2110)
NoiseReductFilterDescTag   = (0x0045,0x103B)
InternalReconAlgorithmTag  = (0x0045,0x1055)

#------------------------
# DICOM hierarchy classes
#------------------------

# DICOM hierarchy levels
DicomPatientLevel     = 0
DicomStudyLevel       = 1
DicomSeriesLevel      = 2
DicomAcquisitionLevel = 3

class Formats(object):
    """
    Formatting options.
    """
    Raw    = 0 # The format specified by DICOM
    Type   = 1 # datetime.date for dates and datetime.time for times
    String = 2 # A string suitable for display
