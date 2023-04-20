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
Class for encapsulating DICOM data sets and extracting DICOM attributes.

"""

import collections

import gdcm

class DicomData(object):
    """
    Class for encapsulating DICOM data sets and extracting DICOM attributes.
    """
    
    def __init__(self, reader=None, file=None, dataSet=None, strFilter=None):
        """
        Optionally initialize the object's attributes.
        """
        
        # Initialize attributes
        self.reader = reader       # The file reader
        self.file = file           # The file object
        self.dataSet = dataSet     # The data set object
        self.strFilter = strFilter # The string filter
         
    def GetAttribute(self, groupElem):
        """
        Return a specified DICOM attribute.

        Parameters
        ----------
        groupElem : tuple, or sequence of tuples
          Either a pair of integers representing the group and element of 
          the desired attribute, or a sequence of such pairs to reach nested 
          attributes. Attributes are converted to a type appropriate for 
          their VR. Sequences are returned as DicomData objects.
        """
        
        # Make sure we've read something
        if self.dataSet is None:
            raise RuntimeError('No DICOM file has been loaded')
         
        # Is the argument a group/element tuple?
        def isGroupElem(ge): return (isinstance(ge, tuple) and len(ge) == 2 
            and isinstance(ge[0], int) and isinstance(ge[1], int))

        # Make sure the input is a group/elem or a list of group/elem tuples
        if not isGroupElem(groupElem) and (not isinstance(groupElem, collections.
            Sequence) or len([ge for ge in groupElem if not isGroupElem(ge)]) != 0):
            raise TypeError("Expected a list of group/elem tuples")

        # Turn single (group, elem) into list
        if isGroupElem(groupElem): groupElem = [groupElem]

        # Handle edge cases gracefully
        if len(groupElem) == 0: return None
        
        # Mapping from VR to Python type function
        vrMap = { "US":int, "UL":int, "SS":int, "SL":int, 
            "FL":float, "FD":float, "DS":float, "IS":int }

        def RecursiveFunc(results, dataSet, groupElem, depth):
            """
            Local recursive worker function
            """

            # Check attribute is present
            ge = groupElem[depth]
            tag = gdcm.Tag(ge[0], ge[1])
            if not dataSet.FindDataElement(tag): return None
 
            # Get the data element
            de = dataSet.GetDataElement(tag)
            
            # Characterize the stage we are at
            isSeq = (str(de.GetVR()) == "SQ")
            atEnd = (depth == len(groupElem) - 1)
            
            # If current attribute is a sequence
            if isSeq:
                
                # Get the sequence of items
                sq = de.GetValueAsSQ()
                nItems = sq.GetNumberOfItems()

                # Loop over each one
                for i in range(nItems):

                    # Get the nested data set
                    ds = sq.GetItem(i + 1).GetNestedDataSet()

                    # Save data set if we at the end, otherwise make the recursive call
                    if atEnd: results.append(DicomData(self.reader, self.file, ds, self.strFilter))
                    else: RecursiveFunc(results, ds, groupElem, depth + 1)
                
            # else if at end of (group, elem) list
            elif atEnd:
                
                # Get the attribute as string
                attrStr = self.strFilter.ToString(de).strip()

                # If there is anything in the string
                if attrStr:
                
                    # Split into component parts
                    attrList = attrStr.split("\\")

                    # Separate and parse to appropriate type
                    vr = str(de.GetVR())
                    if vr in vrMap: attrList = [vrMap[vr](attr) for attr in attrList]
                
                    # Append attribute to results
                    results.append(attrList if len(attrList) > 1 else attrList[0])

            else: raise ValueError("Expected a DICOM sequence to traverse")
        
        # Initialize results to empty list
        results = []

        # Traverse the list and collect the results
        RecursiveFunc(results, self.dataSet, groupElem, 0)

        # Return the results
        if (len(results) == 1 and not isinstance(
            results[0], DicomData)): return results[0]
        elif len(results) == 0: return None
        else: return results
