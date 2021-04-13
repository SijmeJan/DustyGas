#!/bin/bash

# Default: no periodic boundaries, only parameter is exahype file
periodic=false
exahype_file=$1

# Check for flag -p indicating we want periodic boundaries
while getopts ":p" opt; do
  case ${opt} in
    p )
      # Check if mapping file exist, indicating that allow_periodic.sh was run
      FILE=./ExaHyPE-Engine/ExaHyPE/exahype/mappings/PlotPeriodic.cpp
      if [ -f "$FILE" ]; then
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
