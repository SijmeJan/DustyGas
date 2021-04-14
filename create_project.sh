#!/bin/bash

# Exit when error occurs
set -e

# Default: no periodic boundaries, only parameter is exahype file
periodic=false
exahype_file=$1

# Check for flag -p indicating we want periodic boundaries
while getopts ":ph" opt; do
    case ${opt} in
        h )
            echo "Create a dusty gas project from an ExaHyPE file."
            echo
            echo "Usage:"
            echo "    create_project [-p] exahype_file"
            echo
            echo "Options:"
            echo "    -p    Use periodic domain"
            exit 0
            ;;
        p )
            # Check if mapping file exist
            FILE=./ExaHyPE-Engine/ExaHyPE/exahype/mappings/PlotPeriodic.cpp
            if [ -f "$FILE" ]; then
                # Mapping file exists, all is fine
                periodic=true
            else
                # Need to set up periodic boundaries
                ./allow_periodic.sh
                periodic=true
            fi
            echo "Using periodic boundaries" >& 2
            # Second argument should be exahype file
            exahype_file=$2
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: create_project [-d] exahype_file" >&2
            exit 1
            ;;
    esac
done

# Export COMPILER, DISTRIBUTEDMEM
export COMPILER=GNU
export DISTRIBUTEDMEM=MPI

# Use Toolkit to generate source code
echo Generating generic source code...
ExaHyPE-Engine/Toolkit/toolkit.sh $exahype_file

# Write source code
echo Writing specific solver files...
python python/write_cpp --periodic $exahype_file

# Get output directory from exahype file
output_string=$(grep output-directory $exahype_file)
IFS=' ' read -ra ADDR <<< "$output_string"
for i in "${ADDR[@]}"; do
  output_dir=$i
done

output_dir=$(dirname "$exahype_file")"/"$output_dir

echo Starting Make...
cd $output_dir
make
