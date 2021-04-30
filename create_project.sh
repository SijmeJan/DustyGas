#!/bin/bash

# Exit when error occurs
set -e

# Default: no periodic boundaries, only parameter is INI file
periodic=false
ini_file=$1

# Check for flag -p indicating we want periodic boundaries
while getopts ":ph" opt; do
    case ${opt} in
        h )
            echo "Create a dusty gas project from an INI file."
            echo
            echo "Usage:"
            echo "    create_project [-p] INI_FILE"
            echo
            echo "Options:"
            echo "    -p    Use periodic domain"
            exit 0
            ;;
        p )
            # User indicated periodic boundaries; need to check if the
            # necessary modifications have been made. We do this by checking
            # if the PlotPeriodic mapping file exists.
            FILE=./ExaHyPE-Engine/ExaHyPE/exahype/mappings/PlotPeriodic.cpp
            if [ -f "$FILE" ]; then
                # Mapping file exists, all is assumed to be fine
                periodic=true
            else
                # Need to set up periodic boundaries
                ./allow_periodic.sh
                periodic=true
            fi
            echo "Using periodic boundaries" >& 2
            # Second argument should be exahype file
            ini_file=$2
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: create_project [-d] INI_FILE" >&2
            exit 1
            ;;
    esac
done

# Find suitable python on the path
has() { type $@ &>/dev/null; } # a way to check if command is available
if has python3; then PYTHON3="python3";
elif python --version | grep -qi "python 3"; then PYTHON3="python"
else echo "$0: Python3 required for running the ExaHyPE toolkit" >&2; exit -1; fi

# Create Exahype file from ini file
if $periodic
then
    $PYTHON3 python/write_cpp/create_exahype.py --periodic $ini_file
else
    $PYTHON3 python/write_cpp/create_exahype.py $ini_file
fi

# Now work with exahype file: same name, different extension
exahype_file=${ini_file%.ini}".exahype"

# Get output directory from exahype file, where we want to put the executable.
# Note: this directory is relative to the ExaHyPE directory!
output_string=$(grep output-directory $exahype_file)
IFS=' ' read -ra ADDR <<< "$output_string"
for i in "${ADDR[@]}"; do
  output_dir=$i
done
output_dir=$(dirname "$exahype_file")"/"$output_dir

# Copy log setup to output directory
cp $(dirname "$exahype_file")"/exahype.log-filter" $output_dir

# Use GNU compiler
export COMPILER=GNU
#export DISTRIBUTEDMEM=MPI

# Use Toolkit to generate source code
echo Generating generic source code...
ExaHyPE-Engine/Toolkit/toolkit.sh $exahype_file

# Write source code
echo Writing specific solver files...
if $periodic
then
   $PYTHON3 python/write_cpp --periodic $exahype_file
else
   $PYTHON3 python/write_cpp $exahype_file
fi

# Move exahype file into output directory
mv $exahype_file $output_dir"/"

echo Starting Make...
cd $output_dir
make
