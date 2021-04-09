# Restore potentially modified files
#cp backup/ExaHyPE/exahype/runners/Runner.cpp ExaHyPE-Engine/ExaHyPE/exahype/runners/
#cp backup/ExaHyPE/exahype/solvers/Solver.h ExaHyPE-Engine/ExaHyPE/exahype/solvers/
#cp backup/ExaHyPE/exahype/mappings/Prediction.cpp ExaHyPE-Engine/ExaHyPE/exahype/mappings/
#cp backup/ExaHyPE/exahype/solvers/ADERDGSolver.cpp ExaHyPE-Engine/ExaHyPE/exahype/solvers/

cp backup/ExaHyPE/exahype/records/RepositoryState.h ExaHyPE-Engine/ExaHyPE/exahype/records/
cp backup/ExaHyPE/exahype/records/RepositoryState.cpp ExaHyPE-Engine/ExaHyPE/exahype/records/
cp backup/ExaHyPE/exahype/repositories/Repository.h ExaHyPE-Engine/ExaHyPE/exahype/repositories/
cp backup/ExaHyPE/exahype/repositories/RepositorySTDStack.h ExaHyPE-Engine/ExaHyPE/exahype/repositories/
cp backup/ExaHyPE/exahype/repositories/RepositorySTDStack.cpp ExaHyPE-Engine/ExaHyPE/exahype/repositories/
cp backup/ExaHyPE/exahype/runners/Runner.cpp ExaHyPE-Engine/ExaHyPE/exahype/runners/
cp backup/ExaHyPE/exahype/solvers/ADERDGSolver.cpp ExaHyPE-Engine/ExaHyPE/exahype/solvers/
cp backup/ExaHyPE/exahype/solvers/ADERDGSolver.h ExaHyPE-Engine/ExaHyPE/exahype/solvers/
cp backup/ExaHyPE/exahype/solvers/Solver.h ExaHyPE-Engine/ExaHyPE/exahype/solvers/

cp cpp/adapters/* ExaHyPE-Engine/ExaHyPE/exahype/adapters/
cp cpp/mappings/* ExaHyPE-Engine/ExaHyPE/exahype/mappings/

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
make clean
make
