#!/bin/sh

# Directory of IMGT/HLA github repository
wmdadir=$1

# Database version
branch=$1

git --git-dir=$wmdadir/.git --work-tree=$wmdadir checkout $branch

cp $wmdadir/xml/hla.xml.zip .

unzip hla.xml.zip
rm -rf hla.xml.zip

perl validate-annoations.pl hla.xml $branch

