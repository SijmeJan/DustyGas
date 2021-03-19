# Update/generate .exahype file
# Export COMPILER, DISTRIBUTEDMEM
export COMPILER=GNU
export DISTRIBUTEDMEM=MPI

# Use Toolkit to generate source code
echo Generating generic source code...
ExaHyPE-Engine/Toolkit/toolkit.sh $@

# Write source code
echo Modifying fluxes and eigenvalues...
python python/write_cpp $@

# Get output directory from exahype file
output_string=$(grep output-directory $@)
IFS=' ' read -ra ADDR <<< "$output_string"
for i in "${ADDR[@]}"; do
  output_dir=$i
done

output_dir=$(dirname "$@")"/"$output_dir

echo Starting Make...
cd $output_dir
make
