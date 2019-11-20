#!/bin/bash
set -e

target=AlphaPlantImpute2

# Check for tinyhouse
if [[ ! -f alphaplantimpute2/tinyhouse/Pedigree.py ]] ; then
    echo Pedigree.py file not found. Check that the tinyhouse submodule is up to date
    exit 
fi

# Create python wheel
rm -r build dist
python setup.py bdist_wheel

# User guide
cd docs
make latexpdf
cd ..

# Zip file
rm -rf zip
mkdir zip/$target

# Copy the wheel
cp dist/* zip/$target 

# Copy the manual
cp docs/build/latex/*.pdf zip/$target 

# Copy the examples
cp -r example zip/$target 
cd zip
zip -r $target.zip $target
cd ..
