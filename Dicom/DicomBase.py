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
Base and derived classes for the DICOM hierarchy. Trees made
up of these objects are created by DicomScanner, and used by
DicomChooser.
"""

import numpy as np
import datetime


from Dicom.DicomConstants import *

class DicomBase(object):
    """
    Base class for the DICOM hierarchy objects. Holds the attributes
    for the object and defines the GetAttribute function for accessing
    them.
    """
    
    def __init__(self):
        """
        Initialize the base class attributes.
        """
        self.parent = None
        self.children = []
        self.attrs = {}

    def __iter__(self):
        """
        Provide direct iterable access to list of child objects.
        """
        return iter(self.children)

    def GetLevel(self):
        """
        Return the level of this object in the DICOM hierarchy.
        To be overridden by derived classes.
        """
        return None

    def GetParent(self):
        """
        Return the parent object, or None if no parent.
        """
        return self.parent() if self.parent is not None else None

    def GetChildren(self):
        """
        Return the list of child objects.
        """
        return self.children

    def GetAttribute(self, groupElem):
        """
        Return an attribute given its (group,element) tag, or 
        None if it is not present in the attributes dict.
        """
        return self.attrs[groupElem] if groupElem in self.attrs else None

#---------------------------------

class DicomAcquisition(DicomBase):
    """
    Class holding data for a DICOM acquisition. The children of
    the acquisition are the individual instance file names.
          
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
    """

    def __init__(self):
        """
        Initialize the object.
        """
        
        # Call the base class function
        super(DicomAcquisition, self).__init__()

        # Initilize the object's attributes
        self.notes = ""
        self.dimensions = None
        self.spacing = None
        self.homogeneous = None
        self.acquisitionAxis = None
        self.imgPosPats = None
        self.sopInstUIDs = None

    def GetLevel(self):
        """
        Return the level of this object in the DICOM hierarchy.
        """
        return DicomAcquisitionLevel

    def GetFilePaths(self):
        """
        Return the list of child file paths. Functionally identical
        to base class GetChildren.
        """
        return self.children

    def GetNotes(self, format=Formats.Raw):
        """
        Return the user-specified notes.
        """
        if format == Formats.Raw or self.notes is not None: return self.notes
        return ""

    def GetImageType(self, format=Formats.Raw):
        """
        Return the entries of the Image Type attribute, or None if not present.
        """
        attr = self.attrs[ImageTypeTag] if ImageTypeTag in self.attrs else None
        if format == Formats.Raw: return attr
        return ", ".join(attr) if attr is not None else "No Image Types"

    def GetImageComments(self, format=Formats.Raw):
        """
        Return the Image Comments attribute, or None if not present.
        """
        attr = self.attrs[ImageCommentsTag] if ImageCommentsTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return ""

    def GetProtocolName(self, format=Formats.Raw):
        """
        Return the Protocol Name attribute, or None if not present.
        """
        attr = self.attrs[ProtocolNameTag] if ProtocolNameTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return ""

    def GetAcquisitionDate(self, format=Formats.Raw):
        """
        Return the Acquisition Date attribute, or None if not present. If raw
        is False, the date will be converted to a datetime.date and returned.
        """
        attr = self.attrs[AcquisitionDateTag] if AcquisitionDateTag in self.attrs else None
        if format == Formats.Raw: return attr
        if attr is None: return datetime.date if format == Formats.Type else "No Acquis. Date"
        acqDate = datetime.date(int(attr[:4]), int(attr[4:6]), int(attr[6:]))
        if format == Formats.Type: return acqDate
        return acqDate.strftime("%Y-%b-%d")

    def GetAcquisitionTime(self, format=Formats.Raw):
        """
        Return the Acquisition Time attribute, or None if not present. If raw
        is False, the date will be converted to a datetime.time and returned.
        """
        attr = self.attrs[AcquisitionTimeTag] if AcquisitionTimeTag in self.attrs else None
        if format == Formats.Raw: return attr
        if attr is None: return datetime.time if format == Formats.Type else "No Acquis. Time"
        if len(attr) <= 6: acqTime = datetime.time(int(attr[:2]), int(attr[2:4]), int(attr[4:6]))
        else: acqTime = datetime.time(int(attr[:2]), int(attr[2:4]), int(attr[4:6]), int(round(1000.0*float(attr[6:]))))
        if format == Formats.Type: return acqTime
        return acqTime.strftime("%H:%M:%S")

    def GetAcquisitionDateTime(self):
        """
        Return a datetime object constructed by calling GetAcquisitionDate
        and GetAcquisitionTime.
        """
        acq_date = self.GetAcquisitionDate(Formats.Type)
        acq_time =  self.GetAcquisitionTime(Formats.Type)
        return datetime.datetime(acq_date.year, acq_date.month, acq_date.day,
                         acq_time.hour, acq_time.minute, acq_time.second)

    def GetDimensions(self, format=Formats.Raw):
        """
        Return the acquisition dimensions in z,y,x order.
        """
        if format == Formats.Raw: return self.dimensions
        if self.dimensions is None: return "No Dimensions"
        if len(self.dimensions) == 2: return "({} x {})".format(*self.dimensions)
        return "({} x {} x {})".format(*self.dimensions)

    def GetSpacing(self, format=Formats.Raw):
        """
        Return the acquisition spacing in z,y,x order. If raw is False, the 
        spacing will be returned as a string with two digits of precision.
        """
        if format == Formats.Raw: return self.spacing
        if self.spacing is None: return "No Spacing"
        if len(self.spacing) == 2: 
            if self.spacing[0] is None: return "No spacing"
            else: return "({:.2f} x {:.2f})".format(*self.spacing)
        if self.spacing[0] is None: 
            if self.spacing[1] is None: return "(NA x NA x NA)"
            else: return "(NA x {:.2f} x {:.2f})".format(*self.spacing[1:])
        else:
            if self.spacing[1] is None: return "({:.2f} x NA x NA)".format(self.spacing[0])
            else: return "({:.2f} x {:.2f} x {:.2f})".format(*self.spacing)

    def GetOrigin(self, format=Formats.Raw):
        """
        Return the acquisition origin in z,y,x order. If raw is False, the 
        spacing will be returned as a string with two digits of precision.
        """
        if format == Formats.Raw: return self.imgPosPats[0]
        if self.imgPosPats is None: return "No Origin"
        if len(self.imgPosPats[0]) == 2: return "({:.2f} x {:.2f})".format(*self.imgPosPats[0])
        return "({:.2f} x {:.2f} x {:.2f})".format(*self.imgPosPats[0])

    def GetHomogeneous(self, format=Formats.Raw):
        """
        Return a bool indicating whether or not the spacing between slices
        of the acquisition is homogeneous. If the acquisition has only one
        slice, None will be returned.
        """
        if format == Formats.Raw: return self.homogeneous
        if self.homogeneous is None: return "2D"
        return "True" if self.homogeneous else "False"

    def GetAcquisitionAxis(self, format=Formats.Raw):
        """
        Return the acquisition axis in z,y,x order. This is the cross-product
        of the first three elements of Image Orientation (Patient) (0020,0037)
        with the last three elements. If raw is False, the axis will be returned
        as a string with the name (Axial, Sagittal, Coronal) prepended.
        """

        # If a formatted output is desired
        if format != Formats.Raw:

            # If no axis data
            if self.acquisitionAxis is None: return "No Acquis. Axis"

            # Get the truncated string (don't round or vector won't be normalized)
            aa0 = self.acquisitionAxis[0]; as0 = "{:.6f}".format(aa0); d0 = as0.find(".") + 3
            aa1 = self.acquisitionAxis[1]; as1 = "{:.6f}".format(aa1); d1 = as1.find(".") + 3
            aa2 = self.acquisitionAxis[2]; as2 = "{:.6f}".format(aa2); d2 = as2.find(".") + 3
            acqVecStr = "({} x {} x {})".format(as0[:d0], as1[:d1], as2[:d2])

            # Get the dominant axis
            absAcqVec = np.abs(np.array(self.acquisitionAxis))
            iMaxAbsAcqVec = np.argmax(absAcqVec)

            # If it's dominant enough
            if absAcqVec[iMaxAbsAcqVec] > 0.95:

                # Prepend a description and return
                if iMaxAbsAcqVec == 0: acqVecDescStr = "Axial"
                elif iMaxAbsAcqVec == 1: acqVecDescStr = "Coronal"
                else: acqVecDescStr = "Sagittal"
                return "{} ({})".format(acqVecDescStr, self.acquisitionAxis[iMaxAbsAcqVec])

            # Otherwise just return the string
            else: return acqVecStr

        # Just output the raw vector
        return self.acquisitionAxis

    def GetImagePositionPatients(self):
        """
        Return the sorted list of Image Position (Patient) attributes.
        """
        return self.imgPosPats

    def GetSOPInstanceUIDs(self):
        """
        Return the sorted list of SOP Instance UIDs.
        """
        return self.sopInstUIDs

#---------------------------------

class DicomSeries(DicomBase):
    """
    Class holding data for a DICOM series.   
    """

    def __init__(self):
        """
        Initialize the object.
        """
        
        # Call the base class function
        super(DicomSeries, self).__init__()

    def GetLevel(self):
        """
        Return the level of this object in the DICOM hierarchy.
        """
        return DicomSeriesLevel

    def GetAcquisitions(self):
        """
        Return the list of child acquisitions. Functionally identical
        to base class GetChildren.
        """
        return self.children

    def GetSeriesInstanceUID(self, format=Formats.Raw):
        """
        Return the Series Instance UID attribute, or None if not present.
        """
        attr = self.attrs[SeriesInstanceUIDTag] if SeriesInstanceUIDTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return ""

    def GetSeriesDescription(self, format=Formats.Raw):
        """
        Return the Series Description attribute, or None if not present.
        """
        attr = self.attrs[SeriesDescriptionTag] if SeriesDescriptionTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return "No Series Description"

    def GetModality(self, format=Formats.Raw):
        """
        Return the Modality attribute, or None if not present.
        """
        attr = self.attrs[SeriesModalityTag] if SeriesModalityTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return "NA"

    def GetFrameOfReferenceUID(self, format=Formats.Raw):
        """
        Return the Frame Of Reference UID attribute, or None if not present.
        """
        attr = self.attrs[FrameOfReferenceUIDTag] if FrameOfReferenceUIDTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return ""

#---------------------------------

class DicomStudy(DicomBase):
    """
    Class holding data for a DICOM study.   
    """

    def __init__(self):
        """
        Initialize the object.
        """
        
        # Call the base class function
        super(DicomStudy, self).__init__()

    def GetLevel(self):
        """
        Return the level of this object in the DICOM hierarchy.
        """
        return DicomStudyLevel

    def GetSeries(self):
        """
        Return the list of child series. Functionally identical
        to base class GetChildren.
        """
        return self.children

    def GetStudyInstanceUID(self, format=Formats.Raw):
        """
        Return the Study Instance UID attribute, or None if not present.
        """
        attr = self.attrs[StudyInstanceUIDTag] if StudyInstanceUIDTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return ""

    def GetStudyDescription(self, format=Formats.Raw):
        """
        Return the Study Description attribute, or None if not present.
        """
        attr = self.attrs[StudyDescriptionTag] if StudyDescriptionTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return "No Study Description"

    def GetStudyDate(self, format=Formats.Raw):
        """
        Return the Study Date attribute, or None if not present. If raw is
        False, the date will be converted to a datetime.date and returned.
        """
        attr = self.attrs[StudyDateTag] if StudyDateTag in self.attrs else None
        if format == Formats.Raw: return attr
        if attr is None: return "No Study Date"
        acqDate = datetime.date(int(attr[:4]), int(attr[4:6]), int(attr[6:]))
        if format == Formats.Type: return acqDate
        return acqDate.strftime("%Y-%b-%d")

    def GetAccessionNumber(self, format=Formats.Raw):
        """
        Return the Accession Number attribute, or None if not present.
        """
        attr = self.attrs[AccessionNumberTag] if AccessionNumberTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return "No Accession #"

#---------------------------------

class DicomPatient(DicomBase):
    """
    Class holding data for a DICOM patient.   
    """

    def __init__(self):
        """
        Initialize the object.
        """
        
        # Call the base class function
        super(DicomPatient, self).__init__()

    def GetLevel(self):
        """
        Return the level of this object in the DICOM hierarchy.
        """
        return DicomPatientLevel

    def GetStudies(self):
        """
        Return the list of child studies. Functionally identical
        to base class GetChildren.
        """
        return self.children

    def GetPatientsName(self, format=Formats.Raw):
        """
        Return the Patient's Name attribute, or None if not present.
        """
        attr = self.attrs[PatientsNameTag] if PatientsNameTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return "No Patient's Name"

    def GetPatientID(self, format=Formats.Raw):
        """
        Return the Patient ID attribute, or None if not present.
        """
        attr = self.attrs[PatientIDTag] if PatientIDTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return "No Patient ID"

    def GetPatientsAge(self, format=Formats.Raw):
        """
        Return the Patient's Age attribute, or None if not present.
        """
        attr = self.attrs[PatientsAgeTag] if PatientsAgeTag in self.attrs else None
        if format != Formats.Type: return attr if format == Formats.Raw or attr is not None else "NA"
        if attr and len(attr) == 4 and attr[3] in "YMWD":
            try:
                age = float(attr[:3])
                if attr[3] == "D": age /= 365.0
                elif attr[3] == "W": age /= 52.0
                elif attr[3] == "M": age /= 12.0
            except ValueError:
                age = 0.0
        else: age = 0.0
        return age

    def GetPatientsSex(self, format=Formats.Raw):
        """
        Return the Patient's Sex attribute, or None if not present.
        """
        attr = self.attrs[PatientsSexTag] if PatientsSexTag in self.attrs else None
        if format == Formats.Raw or attr is not None: return attr
        return "NA"
