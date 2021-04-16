import os

from common import replace_with_indent, remove_function_body, add_function_body, add_class_member

def allow_periodic(repo_dir):
    ''' Modify ExaHyPE core to allow periodic boundaries on a regular mesh
    '''

    ##########################################################################
    # STEP 1: We need two extra adapters: PlotPeriodic (from mesh to boundary
    # array) and AdjustPeriodic (from boundary array to mesh). The relevant
    # cpp and h files are copied from the cpp directory of the main repo. Here
    # we modify the code to be able to use these adapters.
    ##########################################################################

    # Add adapter names to RepositoryState.h
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/records/RepositoryState.h'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    for i in range(0, len(lines)):
        lines[i] = replace_with_indent(lines[i], 'NumberOfAdapters = 21', 'WriteCheckpoint = 0, ReadCheckpoint = 1, Terminate = 2, RunOnAllNodes = 3, UseAdapterUniformRefinement = 4, UseAdapterMeshRefinement = 5, UseAdapterMeshRefinementAndPlotTree = 6, UseAdapterFinaliseMeshRefinement = 7, UseAdapterFinaliseMeshRefinementOrLocalRollback = 8, UseAdapterInitialPrediction = 9, UseAdapterFusedTimeStep = 10, UseAdapterPredictionRerun = 11, UseAdapterBroadcast = 12, UseAdapterBroadcastAndDropNeighbourMessages = 13, UseAdapterRefinementStatusSpreading = 14, UseAdapterPredictionOrLocalRecomputation = 15, UseAdapterMergeNeighbours = 16, UseAdapterUpdateAndReduce = 17, UseAdapterPrediction = 18, UseAdapterCorrection = 19, UseAdapterAdjustPeriodic = 20, UseAdapterPlotPeriodic = 21, UseAdapterEmpty = 22, NumberOfAdapters = 23\n')

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # Add adapter names to RepositoryState.cpp
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/records/RepositoryState.cpp'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    for i in range(0, len(lines)):
        if (lines[i].find('case UseAdapterCorrection: return "UseAdapterCorrection";') != -1):
                lines[i+1:i+1] = ['      case UseAdapterAdjustPeriodic: return "UseAdapterAdjustPeriodic";\n', '      case UseAdapterPlotPeriodic: return "UseAdapterPlotPeriodic";\n']
        lines[i] = replace_with_indent(lines[i], 'NumberOfAdapters=21', 'return "Action(WriteCheckpoint=0,ReadCheckpoint=1,Terminate=2,RunOnAllNodes=3,UseAdapterUniformRefinement=4,UseAdapterMeshRefinement=5,UseAdapterMeshRefinementAndPlotTree=6,UseAdapterFinaliseMeshRefinement=7,UseAdapterFinaliseMeshRefinementOrLocalRollback=8,UseAdapterInitialPrediction=9,UseAdapterFusedTimeStep=10,UseAdapterPredictionRerun=11,UseAdapterBroadcast=12,UseAdapterBroadcastAndDropNeighbourMessages=13,UseAdapterRefinementStatusSpreading=14,UseAdapterPredictionOrLocalRecomputation=15,UseAdapterMergeNeighbours=16,UseAdapterUpdateAndReduce=17,UseAdapterPrediction=18,UseAdapterCorrection=19,UseAdapterAdjustPeriodic=20,UseAdapterPlotPeriodic=21,UseAdapterEmpty=22,NumberOfAdapters=23)";\n')

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # Add switch functions to Repository.h
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/repositories/Repository.h'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    for i in range(0, len(lines)):
        if (lines[i].find('virtual void switchToCorrection() = 0;') != -1):
            lines[i+1:i+1] = ['    virtual void switchToAdjustPeriodic() = 0;\n',
                              '    virtual void switchToPlotPeriodic() = 0;\n']
        if (lines[i].find('virtual bool isActiveAdapterCorrection() const = 0;') != -1):
            lines[i+1:i+1] = ['    virtual bool isActiveAdapterAdjustPeriodic() const = 0;\n', '    virtual bool isActiveAdapterPlotPeriodic() const = 0;\n']

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # Add functions to RepositorySTDStack.h
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/repositories/RepositorySTDStack.h'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    for i in range(0, len(lines)):
        # Include the new adapter headers
        if (lines[i].find('#include "exahype/adapters/Correction.h"') != -1):
            lines[i+1:i+1] = [' #include "exahype/adapters/AdjustPeriodic.h"\n',
                              ' #include "exahype/adapters/PlotPeriodic.h"\n']
        # Declare new grids
        if (lines[i].find('peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Correction> _gridWithCorrection;') != -1):
            lines[i+1:i+1] = ['    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::AdjustPeriodic> _gridWithAdjustPeriodic;\n', '    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotPeriodic> _gridWithPlotPeriodic;\n']
        # Declare new switch functions
        if (lines[i].find('virtual void switchToCorrection();') != -1):
            lines[i+1:i+1] = ['    virtual void switchToAdjustPeriodic();\n',
                              '    virtual void switchToPlotPeriodic();\n']
        # Declare new IsActive functions
        if (lines[i].find('virtual bool isActiveAdapterCorrection() const;') != -1):
            lines[i+1:i+1] = ['    virtual bool isActiveAdapterAdjustPeriodic() const;\n', '    virtual bool isActiveAdapterPlotPeriodic() const;\n']
        # Declare new time measurement functions
        if (lines[i].find('tarch::timing::Measurement _measureCorrectionCPUTime;') != -1):
            lines[i+1:i+1] = ['    tarch::timing::Measurement _measureAdjustPeriodicCPUTime;\n', '    tarch::timing::Measurement _measurePlotPeriodicCPUTime;\n']
        if (lines[i].find('tarch::timing::Measurement _measureCorrectionCalendarTime;') != -1):
            lines[i+1:i+1] = ['    tarch::timing::Measurement _measureAdjustPeriodicCalendarTime;\n', '    tarch::timing::Measurement _measurePlotPeriodicCalendarTime;\n']


    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # Add functions to RepositorySTDStack.cpp
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/repositories/RepositorySTDStack.cpp'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    for i in range(0, len(lines)):
        if (lines[i].find('_gridWithCorrection(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),') != -1):
            lines[i+1:i+1] = ['  _gridWithAdjustPeriodic(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),\n', '  _gridWithPlotPeriodic(_vertexStack,_cellStack,_geometry,_solverState,domainSize,computationalDomainOffset,_regularGridContainer,_traversalOrderOnTopLevel),\n']
        if (lines[i].find('_gridWithCorrection(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),') != -1):
            lines[i+1:i+1] = ['  _gridWithAdjustPeriodic(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),\n','  _gridWithPlotPeriodic(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),\n']
        if (lines[i].find('_gridWithCorrection.restart') != -1):
            lines[i+1:i+1] = ['  _gridWithAdjustPeriodic.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);\n', '  _gridWithPlotPeriodic.restart(domainSize,domainOffset,domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel);\n']
        if (lines[i].find('_gridWithCorrection.terminate') != -1):
            lines[i+1:i+1] = ['  _gridWithAdjustPeriodic.terminate();\n',
                              '  _gridWithPlotPeriodic.terminate();\n']
        if (lines[i].find('case exahype::records::RepositoryState::UseAdapterCorrection:') != -1):
            lines[i+1:i+1] = ['      case exahype::records::RepositoryState::UseAdapterAdjustPeriodic: watch.startTimer(); _gridWithAdjustPeriodic.iterate(); watch.stopTimer(); _measureAdjustPeriodicCPUTime.setValue( watch.getCPUTime() ); _measureAdjustPeriodicCalendarTime.setValue( watch.getCalendarTime() ); break;\n', '      case exahype::records::RepositoryState::UseAdapterPlotPeriodic: watch.startTimer(); _gridWithPlotPeriodic.iterate(); watch.stopTimer(); _measurePlotPeriodicCPUTime.setValue( watch.getCPUTime() ); _measurePlotPeriodicCalendarTime.setValue( watch.getCalendarTime() ); break;\n']
        if (lines[i].find('void exahype::repositories::RepositorySTDStack::switchToCorrection() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterCorrection); }') != -1):
            lines[i+1:i+1] = [' void exahype::repositories::RepositorySTDStack::switchToAdjustPeriodic() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterAdjustPeriodic); }\n', ' void exahype::repositories::RepositorySTDStack::switchToPlotPeriodic() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotPeriodic); }\n']
        if (lines[i].find('bool exahype::repositories::RepositorySTDStack::isActiveAdapterCorrection() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterCorrection; }') != -1):
            lines[i+1:i+1] = [' bool exahype::repositories::RepositorySTDStack::isActiveAdapterAdjustPeriodic() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterAdjustPeriodic; }\n', ' bool exahype::repositories::RepositorySTDStack::isActiveAdapterPlotPeriodic() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotPeriodic; }\n']

    # Add logging functions
    for i in range(0, len(lines)):
        if (lines[i].find('if (logAllAdapters || _measureCorrectionCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Correction \t |  " << _measureCorrectionCPUTime.getNumberOfMeasurements() << " \t |  " << _measureCorrectionCPUTime.getAccumulatedValue() << " \t |  " << _measureCorrectionCPUTime.getValue()  << " \t |  " << _measureCorrectionCalendarTime.getAccumulatedValue() << " \t |  " << _measureCorrectionCalendarTime.getValue() << " \t |  " << _measureCorrectionCPUTime.toString() << " \t |  " << _measureCorrectionCalendarTime.toString() );') != -1):
            lines[i+1:i+1] = ['    if (logAllAdapters || _measureAdjustPeriodicCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| AdjustPeriodic \t |  " << _measureAdjustPeriodicCPUTime.getNumberOfMeasurements() << " \t |  " << _measureAdjustPeriodicCPUTime.getAccumulatedValue() << " \t |  " << _measureAdjustPeriodicCPUTime.getValue()  << " \t |  " << _measureAdjustPeriodicCalendarTime.getAccumulatedValue() << " \t |  " << _measureAdjustPeriodicCalendarTime.getValue() << " \t |  " << _measureAdjustPeriodicCPUTime.toString() << " \t |  " << _measureAdjustPeriodicCalendarTime.toString() );\n', '    if (logAllAdapters || _measurePlotPeriodicCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotPeriodic \t |  " << _measurePlotPeriodicCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotPeriodicCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotPeriodicCPUTime.getValue()  << " \t |  " << _measurePlotPeriodicCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotPeriodicCalendarTime.getValue() << " \t |  " << _measurePlotPeriodicCPUTime.toString() << " \t |  " << _measurePlotPeriodicCalendarTime.toString() );\n']

    for i in range(0, len(lines)):
        if (lines[i].find('_measureCorrectionCPUTime.erase();') != -1):
            lines[i+1:i+1] = ['   _measureAdjustPeriodicCPUTime.erase();\n',
                              '   _measurePlotPeriodicCPUTime.erase();\n']

    for i in range(0, len(lines)):
        if (lines[i].find('_measureCorrectionCalendarTime.erase();') != -1):
            lines[i+1:i+1] = ['   _measureAdjustPeriodicCalendarTime.erase();\n', '   _measurePlotPeriodicCalendarTime.erase();\n']


    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    ##########################################################################
    # STEP 2: Integrate new adapters into time stepping loop. At the start of
    # a time step, that is, after the correction update: First use the
    # 'plotting' adapter to get the boundary values into and array, followed
    # by the 'adjusting' adapter to get the values on the periodic boundaries
    # on the mesh.
    ##########################################################################

    # Modify time stepping
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/runners/Runner.cpp'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    body = remove_function_body(lines, '::runOneTimeStepWithThreeSeparateAlgorithmicSteps')

    for i in range(0, len(body)):
        if (body[i].find('repository.switchToPrediction(); // Cell onto faces') != -1):
            body[i:i] = ['  repository.switchToPlotPeriodic();\n',
                          '  repository.iterate(1, communicatePeanoVertices);\n',
                          '  repository.switchToAdjustPeriodic();\n',
                          '  repository.iterate(1, communicatePeanoVertices);\n',
                          '\n']

    add_function_body(lines, '::runOneTimeStepWithThreeSeparateAlgorithmicSteps',
                      body)

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    ##########################################################################
    # STEP 3: Declare and define the functions performing the 'plotting' and
    # 'adjusting'. For 'adjusting' we rely on the existing adjustSolution
    # functions. For 'plotting', the final PlotPeriodic function will have to
    # be defined in the user abstract class.
    ##########################################################################

    #################################
    # First deal with ADER-DG solver
    #################################
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/solvers/ADERDGSolver.cpp'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    lines[len(lines):len(lines)] = \
      ['\n',
       'void exahype::solvers::ADERDGSolver::AdjustPeriodic(\n',
       '  const int                                          solverNumber,\n',
       '  CellInfo&                                          cellInfo) {\n',
       '  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);\n',
       '  if ( element != NotFound ) {\n',
       '    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];\n',
       '    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {\n',
       '      double* solverSolution =\n',
       '        static_cast<double*>(cellDescription.getSolution());\n',
       '\n',
       '      const int order = getNodesPerCoordinateAxis() - 1;\n',
       '      const int basisX = order + 1;\n',
       '      const int basisY = order + 1;\n',
       '      const int basisZ = (DIMENSIONS == 3 ? order  : 0 ) + 1;\n',
       '      kernels::index idx_u(basisZ, basisY, basisX, getNumberOfVariables());\n',
       '\n',
       '      // Call the solver-defined AdjustPeriodic function\n',
       '      dfor(i, order + 1) {\n',
       '        AdjustPeriodic(\n',
       '        cellDescription.getOffset(),\n',
       '        cellDescription.getSize(),\n',
       '        i,\n',
       '        solverSolution + idx_u(DIMENSIONS == 3 ? i(2) : 0, i(1), i(0), 0));\n',
       '      }\n',
       '    }\n',
       '  }\n',
       '}\n',
       '\n',
       'void exahype::solvers::ADERDGSolver::PlotPeriodic(\n',
       '  const int                                          solverNumber,\n',
       '  CellInfo&                                          cellInfo) {\n',
       '  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);\n',
       '  if ( element != NotFound ) {\n',
       '    CellDescription& cellDescription = cellInfo._ADERDGCellDescriptions[element];\n',
       '    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {\n',
       '      double* solverSolution =\n',
       '        static_cast<double*>(cellDescription.getSolution());\n',
       '\n',
       '      const int order = getNodesPerCoordinateAxis() - 1;\n',
       '      const int basisX = order + 1;\n',
       '      const int basisY = order + 1;\n',
       '      const int basisZ = (DIMENSIONS == 3 ? order  : 0 ) + 1;\n',
       '      kernels::index idx_u(basisZ, basisY, basisX, getNumberOfVariables());\n',
       '\n',
       '      // Call the solver-defined PlotPeriodic function\n',
       '      dfor(i, order + 1) {\n',
       '        PlotPeriodic(\n',
       '        cellDescription.getOffset(),\n',
       '        cellDescription.getSize(),\n',
       '        i,\n',
       '        solverSolution + idx_u(DIMENSIONS == 3 ? i(2) : 0, i(1), i(0), 0));\n',
       '      }\n',
       '    }\n',
       '  }\n',
       '}\n',
       '\n',
       'void exahype::solvers::ADERDGSolver::FinishPeriodic() {\n',
       '  SendPeriodic();\n',
       '}\n',
       '\n',
       ]

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # Add functions to ADERDGsolver.h
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/solvers/ADERDGSolver.h'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    for i in range(0, len(lines)):
        if (lines[i].find('void updateOrRestrict(') != -1):
            lines[i:i] = ['  void PlotPeriodic(\n',
                          '    const int solverNumber,\n',
                          '    CellInfo& cellInfo) final override;\n',
                          '  void AdjustPeriodic(\n',
                          '    const int solverNumber,\n',
                          '    CellInfo& cellInfo) final override;\n',
                          '  virtual void PlotPeriodic(\n',
                          '    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
                          '    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
                          '    const tarch::la::Vector<DIMENSIONS, int>& pos,\n',
                          '    double* const Q) = 0;\n',
                          '  virtual void AdjustPeriodic(\n',
                          '    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
                          '    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
                          '    const tarch::la::Vector<DIMENSIONS, int>& pos,\n',
                          '    double* Q) = 0;\n',
                          '  void FinishPeriodic() final override;\n',
                          '  virtual void SendPeriodic() = 0;\n',
                          '\n']
            break;

    f = open(fname, "w")
    f.writelines(lines)
    f.close()





    #################################
    # Next: FV solver
    #################################
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/solvers/FiniteVolumesSolver.cpp'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    lines[len(lines):len(lines)] = \
      ['\n',
       'void exahype::solvers::FiniteVolumesSolver::AdjustPeriodic(\n',
       '  const int                                          solverNumber,\n',
       '  CellInfo&                                          cellInfo) {\n',
       '  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);\n',
       '  if ( element != NotFound ) {\n',
       '    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];\n',
       '    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {\n',
       '      double* solverSolution =\n',
       '        static_cast<double*>(cellDescription.getSolution());\n',
       '\n',
       '      kernels::idx3 idx(_nodesPerCoordinateAxis+2*_ghostLayerWidth,\n',
       '                        _nodesPerCoordinateAxis+2*_ghostLayerWidth,\n',
       '                        _numberOfVariables+_numberOfParameters);\n',

       '      for (int i = _ghostLayerWidth; i < _nodesPerCoordinateAxis + _ghostLayerWidth; i++) {\n',
       '        for (int j = _ghostLayerWidth; j < _nodesPerCoordinateAxis + _ghostLayerWidth; j++) {\n',
       '          tarch::la::Vector<DIMENSIONS, int> pos(i, j);\n',
       '          AdjustPeriodic(\n',
       '            cellDescription.getOffset(),\n',
       '            cellDescription.getSize(),\n',
       '            pos,\n',
       '            solverSolution + idx(j, i, 0));\n',
       '        }\n',
       '      }\n',
       '    }\n',
       '  }\n',
       '}\n',
       '\n',
       'void exahype::solvers::FiniteVolumesSolver::PlotPeriodic(\n',
       '  const int                                          solverNumber,\n',
       '  CellInfo&                                          cellInfo) {\n',
       '  const int element = cellInfo.indexOfFiniteVolumesCellDescription(solverNumber);\n',
       '  if ( element != NotFound ) {\n',
       '    CellDescription& cellDescription = cellInfo._FiniteVolumesCellDescriptions[element];\n',
       '    if ( cellDescription.getType()==CellDescription::Type::Leaf ) {\n',
       '      double* solverSolution =\n',
       '        static_cast<double*>(cellDescription.getSolution());\n',
       '\n',
       '      kernels::idx3 idx(_nodesPerCoordinateAxis+2*_ghostLayerWidth,\n',
       '                        _nodesPerCoordinateAxis+2*_ghostLayerWidth,\n',
       '                        _numberOfVariables+_numberOfParameters);\n',

       '      for (int i = _ghostLayerWidth; i < _nodesPerCoordinateAxis + _ghostLayerWidth; i++) {\n',
       '        for (int j = _ghostLayerWidth; j < _nodesPerCoordinateAxis + _ghostLayerWidth; j++) {\n',
       '          tarch::la::Vector<DIMENSIONS, int> pos(i, j);\n',
       '          // Call the solver-defined PlotPeriodic function\n',
       '          PlotPeriodic(\n',
       '            cellDescription.getOffset(),\n',
       '            cellDescription.getSize(),\n',
       '            pos,\n',
       '            solverSolution + idx(j, i, 0));\n',
       '        }\n',
       '      }\n',
       '    }\n',
       '  }\n',
       '}\n',
       '\n',
       'void exahype::solvers::FiniteVolumesSolver::FinishPeriodic() {\n',
       '  SendPeriodic();\n',
       '}\n',
       '\n',
       ]

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # Add functions to FiniteVolumesSolver.h
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/solvers/FiniteVolumesSolver.h'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    for i in range(0, len(lines)):
        if (lines[i].find('void updateOrRestrict(') != -1):
            lines[i:i] = ['  void PlotPeriodic(\n',
                          '    const int solverNumber,\n',
                          '    CellInfo& cellInfo) final override;\n',
                          '  void AdjustPeriodic(\n',
                          '    const int solverNumber,\n',
                          '    CellInfo& cellInfo) final override;\n',
                          '  virtual void PlotPeriodic(\n',
                          '    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
                          '    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
                          '    const tarch::la::Vector<DIMENSIONS, int>& pos,\n',
                          '    double* const Q) = 0;\n',
                          '  virtual void AdjustPeriodic(\n',
                          '    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
                          '    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
                          '    const tarch::la::Vector<DIMENSIONS, int>& pos,\n',
                          '    double* Q) = 0;\n',
                          '  void FinishPeriodic() final override;\n',
                          '  virtual void SendPeriodic() = 0;\n',
                          '\n']
            break;

    f = open(fname, "w")
    f.writelines(lines)
    f.close()



    # Add functions to Solver.h
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/solvers/Solver.h'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    for i in range(0, len(lines)):
        if (lines[i].find('* The nonfused update routine.') != -1):
            lines[i-1:i-1] = \
              ['  virtual void PlotPeriodic(\n',
               '    const int solverNumber,\n',
               '    CellInfo& cellInfo) = 0;\n',
               '  virtual void AdjustPeriodic(\n',
               '    const int solverNumber,\n',
               '    CellInfo& cellInfo) = 0;\n',
               '  virtual void FinishPeriodic() = 0;\n',
               '\n']
            break;

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

if __name__ == "__main__":
    # Full path to cpp files
    repo_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../'

    # execute only if run as a script
    allow_periodic(repo_dir)
