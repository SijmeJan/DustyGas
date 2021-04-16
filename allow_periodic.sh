#!/bin/bash

# By default, enable periodic boundaries
enable_periodic=true

# Check for flag -d, in which case we disable periodic boundaries
while getopts ":dh" opt; do
    case ${opt} in
        h )
            echo "Modify ExaHyPE core to allow periodic boundaries."
            echo
            echo "Usage:"
            echo "    allow_periodic [-d]"
            echo
            echo "Options:"
            echo "    -d    Disable periodic boundaries; switch back to vanilla ExaHyPE"
            exit 0
            ;;
        d )
            echo "Disabling periodic boundaries" >& 2
            enable_periodic=false
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: allow_periodic [-d]" >&2
            exit 1
            ;;
    esac
done

# Restore files from backup
cp backup/ExaHyPE/exahype/records/RepositoryState.h ExaHyPE-Engine/ExaHyPE/exahype/records/
cp backup/ExaHyPE/exahype/records/RepositoryState.cpp ExaHyPE-Engine/ExaHyPE/exahype/records/
cp backup/ExaHyPE/exahype/repositories/Repository.h ExaHyPE-Engine/ExaHyPE/exahype/repositories/
cp backup/ExaHyPE/exahype/repositories/RepositorySTDStack.h ExaHyPE-Engine/ExaHyPE/exahype/repositories/
cp backup/ExaHyPE/exahype/repositories/RepositorySTDStack.cpp ExaHyPE-Engine/ExaHyPE/exahype/repositories/
cp backup/ExaHyPE/exahype/runners/Runner.cpp ExaHyPE-Engine/ExaHyPE/exahype/runners/
cp backup/ExaHyPE/exahype/solvers/ADERDGSolver.cpp ExaHyPE-Engine/ExaHyPE/exahype/solvers/
cp backup/ExaHyPE/exahype/solvers/ADERDGSolver.h ExaHyPE-Engine/ExaHyPE/exahype/solvers/
cp backup/ExaHyPE/exahype/solvers/FiniteVolumesSolver.cpp ExaHyPE-Engine/ExaHyPE/exahype/solvers/
cp backup/ExaHyPE/exahype/solvers/FiniteVolumesSolver.h ExaHyPE-Engine/ExaHyPE/exahype/solvers/
cp backup/ExaHyPE/exahype/solvers/Solver.h ExaHyPE-Engine/ExaHyPE/exahype/solvers/

# Remove additional adapters/mappings
rm -f ExaHyPE-Engine/ExaHyPE/exahype/adapters/PlotPeriodic*
rm -f ExaHyPE-Engine/ExaHyPE/exahype/adapters/AdjustPeriodic*
rm -f ExaHyPE-Engine/ExaHyPE/exahype/mappings/PlotPeriodic*
rm -f ExaHyPE-Engine/ExaHyPE/exahype/mappings/AdjustPeriodic*

# Now ExaHyPE is in its vanilla state. If periodic boundaries are desired,
# modify the ExaHyPE core.
if [ "$enable_periodic" = true ] ; then
    echo Setting up ExaHyPE for periodic boundaries...
    cp cpp/adapters/* ExaHyPE-Engine/ExaHyPE/exahype/adapters/
    cp cpp/mappings/* ExaHyPE-Engine/ExaHyPE/exahype/mappings/

    python python/write_cpp/periodic_exahype.py
fi

# Need a fresh build if modified
echo Make clean
cd ExaHyPE-Engine/ExaHyPE
make clean
