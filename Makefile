# *****************************************************************
#
# IBM Confidential
#
# (C) Copyright IBM Corp. 2019-2021 All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
#
# *****************************************************************
DICOM_VERSION=1.0.3

ARTIFACT_REPOSITORY = https://na.artifactory.swg-devops.com/artifactory/api/pypi/wh-icamb-dev-pypi-local
build_img=wh-imaging-dev-docker-local.artifactory.swg-devops.com/cal-dev-cpu:1.4.0
dockerize=docker run -u $(shell id -u):$(shell id -g) -e HOME=`pwd` --network=host -w ${PWD} -v ${PWD}:${PWD} --rm ${build_img}

all: Dicom/version.py

Dicom/version.py: Makefile
	echo "version=\"${DICOM_VERSION}\"" > $@

check:

version:
	@echo ${DICOM_VERSION}

clean:

dist:
	rm -rf build
	${dockerize} python setup.py bdist_wheel

publish:
	@${dockerize} twine upload --repository-url "${ARTIFACT_REPOSITORY}" --username "${WHI_REPO_USER}" --password ${WHI_REPO_PASSWORD} --disable-progress-bar dist/*.whl

.PHONY: dist
