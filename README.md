# Dicom Tools

This folder contains several classes wrapping the Grassroots DICOM library (GDCM) 
to provide easy, pythonic abstractions for many common DICOM tasks. The following
examples illustrate some of the more common uses.

You will need an installation of GDCM of version 2.6 or later compiled with python
bindings turned on in order to use these classes. You will also need an installation
of PyQt5.

## DicomFile Class

DicomFile defines the methods required to read a DICOM object and query its
attributes.

	from Dicom.DicomFile import DicomFile

    # Load an arbitrary DICOM object from file
    fileName = "1.3.6.1.4.1.12201.119004502355671.1.20160113035634.3.dcm"
    dcmFile = DicomFile(fileName)

	# DicomFile inherits from DicomData, which defines the GetAttribute method. 
	# You pass this method a (group,elem) Dicom tag to indicate which attribute 
	# you want

    # Get the study description 
    studyDesc = dcmFile.GetAttribute((0x0008,0x1030))
    print(studyDesc)
    # Output: CHEST/ABDOMEN/PELVIS

    # GetAttribute() will convert to the right type automatically, and will return
	# None if the attribute doesn't exist. For nested attributes, you specify the 
	# list of DICOM tags leading to the attribute you want

    # Get the bounding box top left hand corner, a nested attribute
    bbTopLeft = dcmFile.GetAttribute([(0x0070,0x0001),(0x0070,0x0008),(0x0070,0x0010)])
    print(bbTopLeft)
    # Output: [302.935, 202.169]

    # An alternate way, asking for attributes one at a time
    sq1 = dcmFile.GetAttribute((0x0070,0x0001))      # For VR = SQ, returns a list of datasets
    sq2 = sq1[0].GetAttribute((0x0070,0x0008))       # For VR = SQ, returns a list of datasets
    bbTopLeft = sq2[0].GetAttribute((0x0070,0x0010))
    print(bbTopLeft)
    # Output: [302.935, 202.169]

	# Note that attributes with a value representations (VR) that defines a sequence 
	# (SQ) may contain an arbitrary number of datasets, so asking for a nested 
	# attribute may actually return several of them

## DicomImage Class

DicomImage inherits from DicomFile, and adds the GetPixelData method to extract the
image data from an image object. The image data is also an attribute, (7FE0,0010), 
but you don't need to know or use this to extract it.

    from Dicom.DicomImage import DicomImage

    # Load the image file
    fileName = "IM000001.dcm"
    dcmImg = DicomImage(fileName)

    # Get the pixel data. The method will automatically apply any rescale slope and 
	# intercept found in the file, and will return a np.float32 numpy array. If you 
	# want the data in its native format (for example to save memory for large objects),
	# pass in native=True
    pixData = dcmImg.GetPixelData()

    # Display it using pyplot
    plt.imshow(pixData, cmap="gray")
    plt.show()

## DicomLoader Class

DicomLoader is an independent class whose main function is to load in sets of DICOM
images from a common series, and assemble them into a 3D volume.

	from Dicom.DicomLoader import DicomLoader

    # Load the series. If you specify a directory name, DicomLoader will scan the 
	# directory for DICOM files. It will take the Series Instance UID from the first
	# file it finds, then load all other images in the directory with the same
	# UID. If you specify a file name, it will read the Series Instance UID from that
	# file, and load all others from the same folder with the same UID. Alternately,
	# you can just pass in a list of files if you already know which ones are in the 
	# series
    seriesPath = "./CT 0.625mm"
    dcmLoad = DicomLoader(seriesPath)

	# DicomLoader performs a number of tests to ensure that the images it is going
	# to load make up a valid three-dimensional object with homogeneous (constant)
	# spacing between each pair of adjacent slices

    # Get the pixel data. As with DicomImage, the function will apply any rescale
	# slope and intercept found in the image files, even if they are different for
	# each image (as with PET scans)
    pixData = dcmLoad.GetPixelData()

	# Get the dimensions, spacing, and origin of the volume (the position of the
	# front, right, bottom-most pixel). These are always returned in Z,Y,X order
	# to match the layout in the numpy array
	shape   = dcmLoad.GetDimensions()
	spacing = dcmLoad.GetSpacing()
	origin  = dcmLoad.GetOrigin()

	# Get attributes from the series. As DicomLoader loads images, it saves one of 
	# them at random so that attributes that are constant throughout the series can
	# be retrieved. Here's an example that extracts the series description
	seriesDesc = dcmLoad.GetImage().GetAttribute((0x0008,0x103E))

	# For several attributes that do vary per-slice, DicomLoader makes the lists of
	# the attributes for each slice available through a set of functions. Lists of
	# the attributes are returned by these functions, sorted into the same order that
	# the corresponding slices appear in the numpy array retruned by GetPixelData
    dcmLoad.GetImagePositionPatients()
    dcmLoad.GetSOPInstanceUIDs()
    dcmLoad.GetRescaleInterceptAndSlopes()
    dcmLoad.GetPixelPaddingValues()
    dcmLoad.GetContrastWindows()

## DicomScan Function

The DicomScan function provides a relatively low-level wrapping of GDCM's very fast
file scanning functionality, to rapidly extract attributes from a large number of 
files. For a higher-level abstraction, see the DicomScanner class further below.
To use DicomScan, you pass in the names of one or more files or directories you would
like to recursively scan, a list of the attributes you are interested in, and a callback
function that will be called once for every DICOM object file found during the scan.

    from Dicom.DicomScan import DicomScan

    # Set the path to scan
    scanPath = "/home/data/Sentara-PS2p0/Standardized/Group01"

    # Dict to hold results
    seriesDict = {}

    # The attributes we want
    serUidTag  = (0x0020,0x000E)
    serDescTag = (0x0008,0x103E)
    imgRowsTag = (0x0028,0x0010)
    imgColsTag = (0x0028,0x0011)
    tags = (serUidTag, serDescTag, imgRowsTag, imgColsTag)

    # The callback function. It is called with the path to the current DICOM file,
	# and a dict mapping the DICOM (group,elem) tags you specified to their values in
	# the current instance, if they were found in this file
    def ProcessFunc(path, tvDict):

        # If the file has all the attributes we are interested in
        if serUidTag in tvDict and serDescTag in tvDict and imgRowsTag in tvDict and imgColsTag in tvDict:

            # Get the series UID
            serUid = tvDict[serUidTag]

            # If we haven't seen this series before
            if serUid not in seriesDict:
       
                # Store the data we want
                seriesDict[serUid] = [1, tvDict[imgRowsTag], tvDict[imgColsTag], path, tvDict[serDescTag]]

            # Otherwise increment the image counter
            else: seriesDict[serUid][0] += 1

    # Execute the scan
    DicomScan(scanPath, tags, ProcessFunc)

This example accumulates a list of the DICOM series it encountered, the DICOM Series 
Description for each one, and the number of files, rows and columns in each one. one
could then display the results as follows.

    # Print out the results
    numTotalImages = 0
    numTotalSeries = 0
    for serData in seriesDict.values():
        print("{}\n    {} x {} x {}".format(serData[4], serData[0], serData[1], serData[2]))
        numTotalImages += serData[0]
        numTotalSeries += 1
    print("{} images in {} series".format(numTotalImages, numTotalSeries))

	# Output:
	Abdomen Without Contrast
		115 x 512 x 512
	Dose Report
		1 x 512 x 512
	Abdomen/Pelvis W/O Contrast
		169 x 512 x 512
	WITH
		446 x 512 x 512
	PE MIP sag
		133 x 512 x 512
	LUNG
		82 x 512 x 512

## DicomScanner Class

The DicomScanner class provides a higher-level wrapping of GDCM's Dicom file scanning
capability. It recursively scans a provided set of files and/or directories and accumulates
requested attributes and data into a patient/study/series/acquisition hierarchy. Note that
this is slightly different from the usual DICOM patient/study/series/image hierarchy.
DicomScanner organizes the images in a series into distinct acquisitions, that share the
same values for key attributes like rows, columns, pixel spacing, orientation in space,
and that form a stack of adjacent slices with the same spacing between them. Sometimes 
multiple acquisitions will be packed into a single series, and DicomLoader will throw an
exception if this occurs, so it is useful to have DicomScanner organize all the files
before trying to load each acqusition.

	# Create the scanner object
	dcmScan = DicomScanner()

	# Add desired attributes. The scanner can be configured to collect different attributes
	# at each level of the hierarchy. Several are included by default, as can be seen in the
	# class __init__ function. In the following example, we add Study ID to be collected at
	# the study level. Since this attribute is usually associated with DICOM study modules, 
	# there is no point collecting it at series or acquisition levels as it will be the same 
	# for all series belonging to a study. there are similar functions for attributes to be
	# collected at the patient, series, and acquisition levels
    dcmScan.AddStudyAttribute((0x0020,0x0010))

	# Now the scan can be performed
	scanPath = "/your/path"
	dcmScan.Scan(scanPath)

	# For scans of large numbers of files, it may take several minutes to perform, so you may
	# wish to save the results for later use. Please use .dsc as the file extension, as this is
	# assumed by the classes decsribed further below
	dcmScan.Save("myScan.dsc")

	# You can load it later using the Load method
	dcmScan.Load("myScan.dsc")

	# The results of the scan are stored internally in a hierarchy of DicomPatient,
	# DicomStudy, DicomSeries, and DicomAcquisition objects, all of which are derived from
	# DicomBase and are defined in the DicomBase module. Each DicomPatient can have multiple
	# DicomStudy objects, which can be accessed vis the GetStudies or GetChildren methods.
	# Each DicomStudy can have multiple DicomSeries objects, which can be accessed via the
	# GetSeries or GetChidren methods. Each DicomSeries can have multiple DicomAcquisition
	# objects, which can be accessed via the GetAcquisitions or GetChildren methods. Each of
	# these objects also has a GetParent method which will return a reference to the object
	# containing it, except for DicomPatient objects which will return None. 
	# 
	# The GetChildren methods are useful when traversing the hierarchy as a graph, so that 
	# the descendants of a node in the hierarchy can be accessed without having to know their
	# type. The GetLevel function of the base class will return for each leve of the hierarchy
	# an integer representing its level, 0 for patients, 1 for studies, 2 for series, and 3 for
	# acquisitions. The DicomAcqusition methods GetFilePaths and GetChildren will return the 
	# list of paths to the DICOM files associated with the acquisition, sorted in order of their
	# spatial position along the acquisition axis. Each object has a set of convenience accessor 
	# functions that can return certain default attributes for their level in a variety of formats. 
	# For attributes that were gathered by the scan but do not have accessor functions (such as
	# the Study ID added in th example above), the GetAttribute function can be used to access
	# the attribute in the same way as described above for the DicomFile class. If the attribute
	# is not present, None will be returned.
	#
	# The easiest way to traverse the hierarchy of DICOM objects is to use the overridden DicomBase
	# __iter__ function. Here is an example that traverse the scan results, locates axial CT 
	# acquisitions with homogeneous spacing and ORIGINAL, PRIMARY image types, loads the DICOM files
	# and accesses the volume's pixel data as a 3D numpy array
	
	# Scan and traverse results
	dcmScan = DicomScanner()
	dcmScan.Scan(scanPath)
	for pat in dcmScan:
		for std in pat:
			for ser in std:
				if ser.GetModality() == "CT":
					for acq in ser:

						# If acquisition has homogeneous slice spacing, has
						# an axial acquisition axis, and has ORIGINAL and
						# PRIMARY Image Types
						if acq.GetHomogeneous() \
							and abs(acq.GetAcquisitionAxis()[0]) > 0.8 \
							and acq.GetImageType()[0] == "ORIGINAL" \
							and acq.GetImageType()[1] == "PRIMARY":

							# Load the acquisition and get pixel data
							try:
								dcmLoad = DicomLoader(acq.GetFilePaths())
								pixelData = dcmLoad.GetPixeldata()
								...
							except Exception:
								# Deal with exception
								...
								continue

	# You may wish to print out a text summary of the hierarchical results
	dcmScan.Print("myResults.txt")

	# Here is an example of the kind of output that is generated
	PATIENT 1 of 27
		Patient ID: 005D01DA7BFD0C0A7483540E52166907
		Patient's Name: 138963
		Patient's Age: 075Y
		Patient's Sex: M
		STUDY 1 of 4
			Study Instance UID: 1.2.840.113654.2.70.1.322915554279783904077766945339697558300
			Study Date: 20160423
			Accession Number: 383
			SERIES 1 of 5
				Series Instance UID: 1.2.840.113654.2.70.1.226116286267300691526809114255586459067
				Series Description: PLAIN
				Modality: CT
				ACQUISITION 1 of 1
					Image Type: ['ORIGINAL', 'PRIMARY', 'AXIAL']
					Dimensions: [470, 512, 512]
					Spacing: [1.0, 0.7949, 0.7949]
					Homogeneous: True
					Acquisition Axis: Axial [1.0, 0.0, 0.0]
			SERIES 2 of 5
				Series Instance UID: 1.2.840.113654.2.70.1.64756047698444127511329991062579494350
				Series Description: IMR1
				Modality: CT
				ACQUISITION 1 of 1
					Image Type: ['ORIGINAL', 'PRIMARY', 'AXIAL']
					Dimensions: [470, 512, 512]
					Spacing: [1.0, 0.7949, 0.7949]
					Homogeneous: True
					Acquisition Axis: Axial [1.0, 0.0, 0.0]

	# Finally, the DicomChooser class described next provides a GUI-based tool for executing
	# scans, saving and loading their results as .dsc files, and exploring the results of a scan
	# in a hierarchical tree structure with various filtering capabilities. The DicomChooser class
	# forms the front end of the DicomViewer class described further below, which allows you to
	# choose multiple scanned acquisitions and load them for viewing

## DicomChooser Class

	# TODO

## DicomViewer Class

	# TODO




