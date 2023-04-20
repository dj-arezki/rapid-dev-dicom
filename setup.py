# *****************************************************************
#
# IBM Confidential
# OCO Source Materials
#
# (C) Copyright IBM Corp. 2021 All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
#
# The source code for this program is not published or otherwise
# divested of its trade secrets, irrespective of what has been
# deposited with the U.S. Copyright Office.
#
# *****************************************************************
import setuptools
import sys

sys.path.append("./")
from Dicom.version import version as dcm_version

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rapid-dev-dicom",
    version=dcm_version,
    author="IBM WHI",
    author_email="pmurugas@in.ibm.com",
    description="Rapid Dev DICOM wrapper",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.ibm.com/WH-Imaging/cal-core",
    packages=["Dicom"],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)
