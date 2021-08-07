#!/bin/bash


# Check if RELEASES is set
if [ -z ${RELEASES+x} ]; then 
    echo "RELEASES is not set. Please specify the release versions to load."; 
else 
    echo "Loading IMGT/HLA release versions: $RELEASES";
fi



