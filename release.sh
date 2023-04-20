#!/bin/bash
# *****************************************************************
#
# IBM Confidential
#
# (C) Copyright IBM Corp. 2021 All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
#
# *****************************************************************
master_branch=main

update_version () {
    pr_branch=$1
    old_version=$2
    new_version=$3

    git push origin --delete $pr_branch
    git branch -D $pr_branch
    git checkout -b $pr_branch
    cat Makefile | sed s/"DICOM_VERSION=$old_version"/"DICOM_VERSION=$new_version"/g > Makefile.tmp
    mv Makefile.tmp Makefile
    git add Makefile
    git commit -m "Update version $old_version -> $new_version"
    git push --set-upstream origin $pr_branch
    git checkout $master_branch
}


git checkout $master_branch && git pull

if [ $? -ne 0 ]
then
   echo Unable to swith to master branch. Please fix local changes
   exit 1
fi

dev_version=`make version`
release_version=`echo $dev_version | sed -E s/.dev.*//g`

grep "rapid-dev-dicom $release_version" NEWS.md

if [ $? -ne 0 ]
then
    echo "No NEWS entry found for rapid-dev-dicom $release_version. Please PR one bofore initiating the release process."
    exit 1
fi

x=`echo $dev_version | cut -d. -f 1`
y=`echo $dev_version | cut -d. -f 2`
z=`echo $dev_version | cut -d. -f 3`

next_dev_version=$x.$y.$(($z+1)).dev0

pr_release_branch=pr_auto_release_`echo $release_version | sed s/'\.'/_/g`
pr_bump_branch=pr_auto_bump_`echo $next_dev_version | sed s/'\.'/_/g`
release_branch=release_"$x"_$y

release_tag=RAPID_DEV_DICOM_`echo $release_version | sed s/'\.'/_/g`

echo current version is: $dev_version
echo creating pr branch for $release_version
update_version $pr_release_branch $dev_version $release_version

echo creating pr branch for $next_dev_version
update_version $pr_bump_branch $dev_version $next_dev_version

echo
echo "================================================================"
echo current version is: $dev_version
echo release version is: $release_version
echo new dev version is: $next_dev_version
echo
echo "please create a first PR for branch $release_branch <- $pr_release_branch"
echo " -> https://github.ibm.com/WH-Imaging/rapid-dev-dicom/compare/$release_branch...$pr_release_branch?expand=1"
echo
echo "please create a second PR for branch $master_branch <- $pr_bump_branch"
echo " -> https://github.ibm.com/WH-Imaging/rapid-dev-dicom/compare/$master_branch...$pr_bump_branch?expand=1"
echo
echo "Do not forget to tag $release_branch after successful build."
echo " git checkout $release_branch"
echo " git pull"
echo " git tag $release_tag"
echo " git push --tags"
echo "================================================================"
