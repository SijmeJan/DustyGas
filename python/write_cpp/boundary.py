from common import replace_with_indent, remove_function_body, add_function_body, add_class_member

def write_outflow_boundary(lines, n_vars):
    # Remove current function body
    remove_function_body(lines, 'boundaryValues')

    periodic = ['  // Outflow boundaries, unimportant if periodic\n']

    for n in range(0, n_vars):
        periodic.extend(['  stateOut[{}] = stateIn[{}];\n'.format(n, n)])
    for n in range(0, n_vars):
        periodic.extend(['  fluxOut[{}] = fluxIn[{}];\n'.format(n, n)])

    add_function_body(lines, 'boundaryValues', periodic)

def write_boundary(lines, n_vars, order, offset, size, solver_name, boundary_name):
    # Strategy: replace function body of mapQuantities to fill the
    # boundaryValues vector.
    remove_function_body(lines, 'mapQuantities')

    x_bound = [offset[0], offset[0] + size[0]]
    y_bound = [offset[1], offset[1] + size[1]]


    boundary = ['  // Fill a boundary array for setting periodic boundaries in 2D, non-AMR runs.\n',
                '  // If mesh = nx times ny, the first nx entries correspond to the bottom boundary.\n',
                '  // The second nx entries correspond to the top boundary.\n',
                '  // The next ny entries correspond to the left boundary.\n',
                '  // The next ny entries correspont to the right boundary.\n',
                '\n',
                '  assertion(outputQuantities==nullptr);\n\n',
                '  // Hack: number of cells in x and y\n',
                '  int n_cell_x = (int) round({}/sizeOfPatch[0]);\n'.format(size[0]),
                '  int n_cell_y = (int) round({}/sizeOfPatch[1]);\n'.format(size[1]),
                '\n',
                '  // Global (uniform) cell resolution\n',
                '  global_dx[0] = sizeOfPatch[0];\n',
                '  global_dx[1] = sizeOfPatch[1];\n',
                '\n',
                '  // Global mesh size\n',
                '  global_n[0] = n_cell_x;\n',
                '  global_n[1] = n_cell_y;\n',
                '\n',
                '  // Make sure vector is of the correct size, and set elements to zero.\n',
                '  // Note that this should happen only the first time this function is called.\n',
                '  int n_bound_cells = 2*(global_n[0] + global_n[1]);\n',
                '  int n_send_per_cell = {};\n'.format(n_vars*(order+1)*(order+1)),
                '  if (boundaryValues_local.size() != n_send_per_cell*n_bound_cells) {\n',
                '    boundaryValues_local.resize(n_send_per_cell*n_bound_cells);\n',
                '    std::fill(boundaryValues_local.begin(), boundaryValues_local.end(), 0.0);\n',
                '  }\n',
                '\n',
                '  // Number of cells to left and bottom\n',
                '  int cell_x = (int) ((offsetOfPatch[0] + 0.5*sizeOfPatch[0])/sizeOfPatch[0]);\n',
                '  int cell_y = (int) ((offsetOfPatch[1] + 0.5*sizeOfPatch[1])/sizeOfPatch[1]);\n',
                '\n',
                '  int indx = -1;\n',
                '\n',
                '  if (cell_y == 1) {\n',
                '    indx = cell_x;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '    //std::cout << "Mapping quantities at x = " << x[0] << ", y = " << x[1] << ", position " << pos[0] << " " << pos[1] << ", index = " << indx << ", arr_index = " << arr_index << ", cell index " << cell_x << " " << cell_y << std::endl;\n',
                '  }\n'
                '  if (cell_y == global_n[1] - 2) {\n',
                '    indx = n_cell_x + cell_x;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '    //std::cout << "Mapping quantities at x = " << x[0] << ", y = " << x[1] << ", position " << pos[0] << " " << pos[1] << ", index = " << indx << ", arr_index = " << arr_index << ", cell index " << cell_x << " " << cell_y << std::endl;\n',
                '  }\n'
                '  if (cell_x == 1) {\n',
                '    indx = 2*n_cell_x + cell_y;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '    //std::cout << "Mapping quantities at x = " << x[0] << ", y = " << x[1] << ", position " << pos[0] << " " << pos[1] << ", index = " << indx << ", arr_index = " << arr_index << ", cell index " << cell_x << " " << cell_y << std::endl;\n',
                '  }\n'
                '  if (cell_x == global_n[0] - 2) {\n',
                '    indx = 2*n_cell_x + n_cell_y + cell_y;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '    //std::cout << "Mapping quantities at x = " << x[0] << ", y = " << x[1] << ", position " << pos[0] << " " << pos[1] << ", index = " << indx << ", arr_index = " << arr_index << ", cell index " << cell_x << " " << cell_y << std::endl;\n',
                '  }\n'
                '\n',
                ]

    add_function_body(lines, 'mapQuantities', boundary)

    # Include solver header
    solver = ['#include "tarch/parallel/Node.h"\n',
              '#include "tarch/parallel/NodePool.h"\n',
              '#include "{}.h"\n'.format(solver_name)]

    # Put it after first include command
    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i+1:i+1] = solver
            break;

    # Init pointers
    solver = ['  boundaryValues = &(solver.periodicBoundaryValues);\n',
              '  global_dx = &(solver.global_dx[0]);\n',
              '  global_n = &(solver.global_n[0]);\n']

    # Put it as first line of first function
    remove_function_body(lines, '::' + boundary_name)
    add_function_body(lines, '::' + boundary_name, solver)

    solver = ['  std::fill(boundaryValues_local.begin(), boundaryValues_local.end(), 0.0);\n',
              '  std::cout << "Starting boundary plotting at t = " << time << std::endl;\n']

    remove_function_body(lines, 'startPlotting')
    add_function_body(lines, 'startPlotting', solver)

    solver = ['#ifdef Parallel\n',
              '  // Should only happen the first time of call\n',
              '  if ((*boundaryValues).size() == 0) {\n',
              '    if (tarch::parallel::Node::getInstance().getNumberOfNodes() > 1) {\n',
              '      // Main rank does not know the boundary size\n',
              '      int s_loc = boundaryValues_local.size();\n',
              '      int s_max = 0;\n',
              '      MPI_Reduce(&s_loc, &s_max, 1, MPI_INT, MPI_MAX, 0,\n',
              '                 tarch::parallel::Node::getInstance().getCommunicator());\n',
              '      // Set correct size so we can use Allreduce\n',
              '      if (tarch::parallel::Node::getInstance().isGlobalMaster()) {\n',
              '        boundaryValues_local.resize(s_max);\n',
              '        std::fill(boundaryValues_local.begin(), boundaryValues_local.end(), 0.0);\n',
              '        (*boundaryValues).resize(s_max);\n',
              '      } else {\n',
              '        (*boundaryValues).resize(boundaryValues_local.size());\n',
              '      }\n',
              '    } else {\n',
              '      (*boundaryValues).resize(boundaryValues_local.size());\n',
              '    }\n',
              '  }\n',
              '  //std::cout << "rank = " << tarch::parallel::Node::getInstance().getRank() << ", sizes: " << boundaryValues_local.size() << " " << (*boundaryValues).size() << std::endl;\n',
              '  std::cout << "Allreduce with rank " << tarch::parallel::Node::getInstance().getRank() << " and communicator " << tarch::parallel::Node::getInstance().getCommunicator() << std::endl;\n',
              '  MPI_Allreduce(&(boundaryValues_local[0]),\n',
              '                &((*boundaryValues)[0]),\n',
              '                boundaryValues_local.size(),\n',
              '                MPI_DOUBLE,\n',
              '                MPI_SUM,\n',
              '                tarch::parallel::Node::getInstance().getCommunicator());\n',
              '#else\n',
              '  (*boundaryValues) = boundaryValues_local;\n',
              '#endif\n']

    remove_function_body(lines, 'finishPlotting')
    add_function_body(lines, 'finishPlotting', solver)


def write_boundary_h(lines):
    # Declare pointers in boundary header file
    boundary = [' private:\n',
                '  std::vector<double> *boundaryValues;\n',
                '  std::vector<double> boundaryValues_local;\n',
                '  double *global_dx;\n',
                '  int *global_n;\n']

    # Do nothing if line already present
    for i in range(0, len(lines)):
        if (lines[i].find(boundary[1]) != -1):
            return

    for i in range(0, len(lines)):
        if (lines[i].find('public:') != -1):
            lines[i:i] = boundary
            break;

    # Make sure to include vector
    boundary = ['#include <vector>\n']

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i:i] = boundary
            break;

def write_solver_h(lines):
    # Declare arrays
    solver = ['    std::vector<double> periodicBoundaryValues;\n',
              '    double global_dx[2];\n',
              '    int global_n[2];\n',
              '    int global_dof_index;\n']

    # Do nothing if line already present
    for i in range(0, len(lines)):
        if (lines[i].find(solver[0]) != -1):
            return

    # Put it just before end of class
    for i in range(0, len(lines)):
        if (lines[i].find('};') != -1):
            lines[i:i] = solver
            break;

    # Make sure to include vector header
    solver = ['#include <vector>\n']

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i:i] = solver
            break;

def write_solver_set_periodic(lines, n_vars, order):
    body = remove_function_body(lines, 'adjustPointSolution')

    periodic = ['  if (t < 0.0) {\n',
                '    // Keep a counter of dofs visited\n',
                '    int max_global_dof_index = global_n[0]*global_n[1]*{}*{};\n'.format(order + 1, order + 1),
                '    if (global_dof_index >= max_global_dof_index)\n',
                '      global_dof_index -= max_global_dof_index;\n',
                '\n',
                '    // Global cell number in x and y direction\n',
                '    int i = (int) (x[0]/global_dx[0]);\n',
                '    int j = (int) (x[1]/global_dx[1]);\n',
                '\n',
                '    int indx = -1;\n',
                '    if (j == 0) indx = global_n[0] + i;\n',
                '    if (j == global_n[1] - 1) indx = i;\n',
                '    if (i == 0) indx = 2*global_n[0] + global_n[1] + j;\n',
                '    if (i == global_n[0] - 1) indx = 2*global_n[0] + j;\n',
                '\n',
                '    if (indx >= 0) {\n',
                '      int dof = global_dof_index % {};\n'.format((order + 1)*(order + 1)),
                '      int arr_index = indx*{} + dof*{};\n'.format(n_vars*(order+1)*(order+1), n_vars),
                '      //if (i == 0 && j == 0) std::cout << "Adjusting point solution at x = " << x[0] << ", y = " << x[1] << ", cell index " << i << " " << j << ", index = " << indx << ", arr_index = " << arr_index << " at t = " << t << std::endl;\n',
                '\n',

                '      for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '        Q[n] = periodicBoundaryValues[arr_index + n];\n',
                '    }\n',
                '    global_dof_index++;\n',

                '  }\n']

    body.extend(periodic)

    add_function_body(lines, 'adjustPointSolution', body)

    remove_function_body(lines, 'init')

    body = ['  global_dof_index = 0;\n']
    add_function_body(lines, 'init', body)

    # Remove current function body
    #remove_function_body(lines, 'boundaryValues')

    #periodic = ['  // Outflow boundaries, unimportant since periodic anyway\n']

    #for n in range(0, n_vars):
    #    periodic.extend(['  stateOut[{}] = stateIn[{}];\n'.format(n, n)])
    #for n in range(0, n_vars):
    #    periodic.extend(['  fluxOut[{}] = fluxIn[{}];\n'.format(n, n)])

    #add_function_body(lines, 'boundaryValues', periodic)

    return

    # Remove current function body
    remove_function_body(lines, 'boundaryValues')

    periodic = ['  // Global cell number in x and y direction\n',
                '  int i = (int) round(x[0]/global_dx[0] - 0.5);\n',
                '  int j = (int) round(x[1]/global_dx[1] - 0.5);\n',
                '\n',
                '  int indx = -1;\n',
                '  // Bottom boundary\n',
                '  if (faceIndex == 2) {\n',
                '    indx = global_n[0] + i;\n',
                '    indx += std::abs(j + 1)*2*(global_n[0] + global_n[1]);\n',
                '  }\n',
                '  // Top boundary\n',
                '  if (faceIndex == 3) {\n',
                '    indx = i;\n',
                '    indx += std::abs(j - global_n[1])*2*(global_n[0] + global_n[1]);\n',
                '  }\n',
                '  // Left boundary\n',
                '  if (faceIndex == 0) {\n',
                '    indx = 2*global_n[0] + global_n[1] + j;\n',
                '    indx += std::abs(i + 1)*2*(global_n[0] + global_n[1]);\n',
                '  }\n',
                '  // Right boundary\n',
                '  if (faceIndex == 1) {\n',
                '    indx = 2*global_n[0] + j;\n',
                '    indx += std::abs(i - global_n[0])*2*(global_n[0] + global_n[1]);\n',
                '  }\n',
                '\n']
    #for i in range(0, n_vars):
    #    periodic.append('  stateOutside[{}] = periodicBoundaryValues[{}*indx + {}];\n'.format(i, n_vars, i))

    #periodic.append('  std::cout << "Setting boundary at x = " << x[0] << ", y = " << x[1] << ", i = " << i << ", j = " << j << " " << indx << std::endl;\n')

    add_function_body(lines, 'boundaryValues', periodic)



##########################
# FINITE VOLUME VERSIONS #
##########################

def write_boundary_fv(lines, n_vars, patch_size, offset, size, solver_name, boundary_name):
    # Strategy: replace function body of mapQuantities to fill the
    # boundaryValues vector.
    remove_function_body(lines, 'mapQuantities')

    x_bound = [offset[0], offset[0] + size[0]]
    y_bound = [offset[1], offset[1] + size[1]]

    boundary = ['  // Fill a boundary array for setting periodic boundaries in 2D, non-AMR runs.\n',
                '  // If FV mesh = nx times ny, the first nx entries correspond to the bottom boundary.\n',
                '  // The second nx entries correspond to the top boundary.\n',
                '  // The next ny entries correspond to the left boundary.\n',
                '  // The next ny entries correspont to the right boundary.\n',
                '\n',
                '  assertion(outputQuantities==nullptr);\n\n',
                '  // Hack: number of patches in x and y\n',
                '  int n_patch_x = (int) round({}/sizeOfPatch[0]);\n'.format(size[0]),
                '  int n_patch_y = (int) round({}/sizeOfPatch[1]);\n'.format(size[1]),
                '\n',
                '  // Global (uniform) resolution\n',
                '  global_dx[0] = sizeOfPatch[0]/{};\n'.format(patch_size),
                '  global_dx[1] = sizeOfPatch[1]/{};\n'.format(patch_size),
                '\n',
                '  // Global mesh size\n',
                '  global_n[0] = n_patch_x*{};\n'.format(patch_size),
                '  global_n[1] = n_patch_y*{};\n'.format(patch_size),
                '\n',
                '  // Make sure vecor is of the correct size, and set elements to zero.\n',
                '  // Note that this should happen only the first time this function is called.\n',
                '  int n_bound_cells = 2*nGhostCells*(global_n[0] + global_n[1]);\n',
                '  if (boundaryValues_local.size() != {}*n_bound_cells) {{\n'.format(n_vars),
                '    boundaryValues_local.resize({}*n_bound_cells);\n'.format(n_vars),
                '    std::fill(boundaryValues_local.begin(), boundaryValues_local.end(), 0.0);\n',
                '  }\n',
                '\n',
                '  // Number of patches to left and bottom\n',
                '  int patch_x = (int) ((offsetOfPatch[0] + 0.5*sizeOfPatch[0])/sizeOfPatch[0]);\n',
                '  int patch_y = (int) ((offsetOfPatch[1] + 0.5*sizeOfPatch[1])/sizeOfPatch[1]);\n',
                '\n',
                '  for (int ng = 0; ng < nGhostCells; ng++) {\n',
                '    int indx = -1;\n',
                '\n',
                '    if (x[1] - 0.5*sizeOfPatch[1] < {} && pos[1] == ng) {{\n'.format(y_bound[0]),
                '      // Bottom boundary\n',
                '      indx = patch_x*{} + pos[0];\n'.format(patch_size),
                '      indx = indx + 2*ng*(global_n[0] + global_n[1]);\n',
                '      int arr_index = {}*indx;\n'.format(n_vars),
                '      for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '        boundaryValues_local[arr_index + n] = Q[n];\n',
                '    }\n',
                '    if (x[1] + 0.5*sizeOfPatch[1] > {} && pos[1] == {} - ng) {{\n'.format(y_bound[1], patch_size-1),
                '      // Top boundary\n',
                '      indx = (n_patch_x + patch_x)*{} + pos[0];\n'.format(patch_size),
                '      indx = indx + 2*ng*(global_n[0] + global_n[1]);\n',
                '      int arr_index = {}*indx;\n'.format(n_vars),
                '      for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '        boundaryValues_local[arr_index + n] = Q[n];\n',
                '    }\n',
                '    if (x[0] - 0.5*sizeOfPatch[0] < {} && pos[0] == ng) {{\n'.format(x_bound[0]),
                '      // Left boundary\n',
                '      indx = (2*n_patch_x + patch_y)*{} + pos[1];\n'.format(patch_size),
                '      indx = indx + 2*ng*(global_n[0] + global_n[1]);\n',
                '      int arr_index = {}*indx;\n'.format(n_vars),
                '      for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '        boundaryValues_local[arr_index + n] = Q[n];\n',
                '    }\n',
                '    if (x[0] + 0.5*sizeOfPatch[0] > {} && pos[0] == {} - ng) {{\n'.format(x_bound[1], patch_size-1),
                '      // Right boundary\n',
                '      indx = (2*n_patch_x + n_patch_y + patch_y)*{} + pos[1];\n'.format(patch_size),
                '      indx = indx + 2*ng*(global_n[0] + global_n[1]);\n',
                '      int arr_index = {}*indx;\n'.format(n_vars),
                '      for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '        boundaryValues_local[arr_index + n] = Q[n];\n',
                '    }\n',
                '  }\n']

    add_function_body(lines, 'mapQuantities', boundary)

    # Include solver header
    solver = ['#include "tarch/parallel/Node.h"\n',
              '#include "tarch/parallel/NodePool.h"\n',
              '#include "{}.h"\n'.format(solver_name)]

    # Put it after first include command
    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i+1:i+1] = solver
            break;

    # Init pointers
    solver = ['  boundaryValues = &(solver.periodicBoundaryValues);\n',
              '  global_dx = &(solver.global_dx[0]);\n',
              '  global_n = &(solver.global_n[0]);\n',
              '  nGhostCells = solver.GhostLayerWidth;\n']

    # Put it as first line of first function
    remove_function_body(lines, '::' + boundary_name)
    add_function_body(lines, '::' + boundary_name, solver)

    solver = ['  std::fill(boundaryValues_local.begin(), boundaryValues_local.end(), 0.0);\n']

    remove_function_body(lines, 'startPlotting')
    add_function_body(lines, 'startPlotting', solver)

    solver = ['#ifdef Parallel\n',
              '  // Should only happen the first time of call\n',
              '  if ((*boundaryValues).size() == 0) {\n',
              '    if (tarch::parallel::Node::getInstance().getNumberOfNodes() > 1) {\n',
              '      // Main rank does not know the boundary size\n',
              '      int s_loc = boundaryValues_local.size();\n',
              '      int s_max = 0;\n',
              '      MPI_Reduce(&s_loc, &s_max, 1, MPI_INT, MPI_MAX, 0,\n',
              '                 tarch::parallel::Node::getInstance().getCommunicator());\n',
              '      // Set correct size so we can use Allreduce\n',
              '      if (tarch::parallel::Node::getInstance().isGlobalMaster()) {\n',
              '        boundaryValues_local.resize(s_max);\n',
              '        std::fill(boundaryValues_local.begin(), boundaryValues_local.end(), 0.0);\n',
              '        (*boundaryValues).resize(s_max);\n',
              '      } else {\n',
              '        (*boundaryValues).resize(boundaryValues_local.size());\n',
              '      }\n',
              '    } else {\n',
              '      (*boundaryValues).resize(boundaryValues_local.size());\n',
              '    }\n',
              '  }\n',
              '  std::cout << "rank = " << tarch::parallel::Node::getInstance().getRank() << ", sizes: " << boundaryValues_local.size() << " " << (*boundaryValues).size() << std::endl;\n',
              '  std::cout << "Allreduce with rank " << tarch::parallel::Node::getInstance().getRank() << " and communicator " << tarch::parallel::Node::getInstance().getCommunicator() << std::endl;\n',
              '  MPI_Allreduce(&(boundaryValues_local[0]),\n',
              '                &((*boundaryValues)[0]),\n',
              '                boundaryValues_local.size(),\n',
              '                MPI_DOUBLE,\n',
              '                MPI_SUM,\n',
              '                tarch::parallel::Node::getInstance().getCommunicator());\n',
              '#else\n',
              '  (*boundaryValues) = boundaryValues_local;\n',
              '#endif\n']

    remove_function_body(lines, 'finishPlotting')
    add_function_body(lines, 'finishPlotting', solver)


def write_boundary_h(lines):
    # Declare pointers in boundary header file
    boundary = [' private:\n',
                '  std::vector<double> *boundaryValues;\n',
                '  std::vector<double> boundaryValues_local;\n',
                '  double *global_dx;\n',
                '  int *global_n;\n',
                '  int nGhostCells;\n']

    # Do nothing if line already present
    for i in range(0, len(lines)):
        if (lines[i].find(boundary[1]) != -1):
            return

    for i in range(0, len(lines)):
        if (lines[i].find('public:') != -1):
            lines[i:i] = boundary
            break;

    # Make sure to include vector
    boundary = ['#include <vector>\n']

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i:i] = boundary
            break;

def write_solver_h_fv(lines):
    # Declare arrays
    solver = ['    std::vector<double> periodicBoundaryValues;\n',
              '    double global_dx[2];\n',
              '    int global_n[2];\n']

    # Do nothing if line already present
    for i in range(0, len(lines)):
        if (lines[i].find(solver[0]) != -1):
            return

    # Put it just before end of class
    for i in range(0, len(lines)):
        if (lines[i].find('};') != -1):
            lines[i:i] = solver
            break;

    # Make sure to include vector header
    solver = ['#include <vector>\n']

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i:i] = solver
            break;

def write_solver_set_periodic_fv(lines, n_vars):
    # Remove current function body
    remove_function_body(lines, 'boundaryValues')

    periodic = ['  // Global cell number in x and y direction\n',
                '  int i = (int) round(x[0]/global_dx[0] - 0.5);\n',
                '  int j = (int) round(x[1]/global_dx[1] - 0.5);\n',
                '\n',
                '  int indx = -1;\n',
                '  // Bottom boundary\n',
                '  if (faceIndex == 2) {\n',
                '    indx = global_n[0] + i;\n',
                '    indx += std::abs(j + 1)*2*(global_n[0] + global_n[1]);\n',
                '  }\n',
                '  // Top boundary\n',
                '  if (faceIndex == 3) {\n',
                '    indx = i;\n',
                '    indx += std::abs(j - global_n[1])*2*(global_n[0] + global_n[1]);\n',
                '  }\n',
                '  // Left boundary\n',
                '  if (faceIndex == 0) {\n',
                '    indx = 2*global_n[0] + global_n[1] + j;\n',
                '    indx += std::abs(i + 1)*2*(global_n[0] + global_n[1]);\n',
                '  }\n',
                '  // Right boundary\n',
                '  if (faceIndex == 1) {\n',
                '    indx = 2*global_n[0] + j;\n',
                '    indx += std::abs(i - global_n[0])*2*(global_n[0] + global_n[1]);\n',
                '  }\n',
                '\n']
    for i in range(0, n_vars):
        periodic.append('  stateOutside[{}] = periodicBoundaryValues[{}*indx + {}];\n'.format(i, n_vars, i))

    #periodic.append('  std::cout << "Setting boundary at x = " << x[0] << ", y = " << x[1] << ", i = " << i << ", j = " << j << " " << indx << std::endl;\n')

    add_function_body(lines, 'boundaryValues', periodic)


def correction_boundary_hack(repo_dir, offset, size, n_vars, order):
    # Add adapter names to RepositoryState
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
        if (lines[i].find('#include "exahype/adapters/Correction.h"') != -1):
            lines[i+1:i+1] = [' #include "exahype/adapters/AdjustPeriodic.h"\n',
                              ' #include "exahype/adapters/PlotPeriodic.h"\n']
        if (lines[i].find('peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Correction> _gridWithCorrection;') != -1):
            lines[i+1:i+1] = ['    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::AdjustPeriodic> _gridWithAdjustPeriodic;\n', '    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotPeriodic> _gridWithPlotPeriodic;\n']
        if (lines[i].find('virtual void switchToCorrection();') != -1):
            lines[i+1:i+1] = ['    virtual void switchToAdjustPeriodic();\n',
                              '    virtual void switchToPlotPeriodic();\n']
        if (lines[i].find('virtual bool isActiveAdapterCorrection() const;') != -1):
            lines[i+1:i+1] = ['    virtual bool isActiveAdapterAdjustPeriodic() const;\n', '    virtual bool isActiveAdapterPlotPeriodic() const;\n']
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

    # Add function to ADERDGsolver.cpp
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
       '    if ( cellDescription.getType()==CellDescription::Type::Leaf )\n',
       '      adjustSolutionAfterUpdate(cellDescription);\n',
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
       '      const int order = getNodesPerCoordinateAxis();\n',
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
       '}\n\n']

    f = open(fname, "w")
    f.writelines(lines)
    f.close()




    # Add function to ADERDGsolver.h
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
                          '    double* const Q) = 0;\n\n']
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
               #'  virtual void PlotPeriodic(\n',
               #'    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
               #'    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
               #'    const tarch::la::Vector<DIMENSIONS, int>&    pos,\n',
               #'    double* const Q) = 0;\n',
               '\n']
            break;

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    return

    # Add function to DustyGasSolver.cpp
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/solvers/ADERDGSolver.cpp'
    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    x_bound = [offset[0], offset[0] + size[0]]
    y_bound = [offset[1], offset[1] + size[1]]

    # TODO: declare periodic array in header
    # TODO: should this be in DustyGasSolver?
    boundary = ['void exahype::solvers::ADERDGSolver::PlotPeriodic(\n',
                '  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
                '  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
                '  const tarch::la::Vector<DIMENSIONS, int>&    pos,\n',
                '  double* const Q\n',
                '  ) {\n',
                '/*',
                '  // Fill a boundary array for setting periodic boundaries in 2D, non-AMR runs.\n',
                '  // If mesh = nx times ny, the first nx entries correspond to the bottom boundary.\n',
                '  // The second nx entries correspond to the top boundary.\n',
                '  // The next ny entries correspond to the left boundary.\n',
                '  // The next ny entries correspont to the right boundary.\n',
                '\n',
                '  // Hack: number of cells in x and y\n',
                '  int n_cell_x = (int) round({}/sizeOfPatch[0]);\n'.format(size[0]),
                '  int n_cell_y = (int) round({}/sizeOfPatch[1]);\n'.format(size[1]),
                '\n',
                '  // Global (uniform) cell resolution\n',
                '  global_dx[0] = sizeOfPatch[0];\n',
                '  global_dx[1] = sizeOfPatch[1];\n',
                '\n',
                '  // Global mesh size\n',
                '  global_n[0] = n_cell_x;\n',
                '  global_n[1] = n_cell_y;\n',
                '\n',
                '  // Make sure vector is of the correct size, and set elements to zero.\n',
                '  // Note that this should happen only the first time this function is called.\n',
                '  int n_bound_cells = 2*(global_n[0] + global_n[1]);\n',
                '  int n_send_per_cell = {};\n'.format(n_vars*(order+1)*(order+1)),
                '  if (boundaryValues_local.size() != n_send_per_cell*n_bound_cells) {\n',
                '    boundaryValues_local.resize(n_send_per_cell*n_bound_cells);\n',
                '    std::fill(boundaryValues_local.begin(), boundaryValues_local.end(), 0.0);\n',
                '  }\n',
                '\n',
                '  // Number of cells to left and bottom\n',
                '  int cell_x = (int) ((offsetOfPatch[0] + 0.5*sizeOfPatch[0])/sizeOfPatch[0]);\n',
                '  int cell_y = (int) ((offsetOfPatch[1] + 0.5*sizeOfPatch[1])/sizeOfPatch[1]);\n',
                '\n',
                '  int indx = -1;\n',
                '\n',
                '  if (cell_y == 1) {\n',
                '    indx = cell_x;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '  }\n'
                '  if (cell_y == global_n[1] - 2) {\n',
                '    indx = n_cell_x + cell_x;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '  }\n'
                '  if (cell_x == 1) {\n',
                '    indx = 2*n_cell_x + cell_y;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '  }\n'
                '  if (cell_x == global_n[0] - 2) {\n',
                '    indx = 2*n_cell_x + n_cell_y + cell_y;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '  }\n'
                '*/',
                '}\n',
                '\n',
                ]

    lines[len(lines):len(lines)] = boundary

    # ExaHyPE fuses the correction step with the 'adjustSolution' step, and the
    # prediction step with the plotting step. For periodic boundaries with
    # ADERDG, this is unwanted. One of the plotters 'plots' the periodic data
    # into an array, while the 'adjustSolution' gets the data from the array
    # to the grid.
    #
    # We therefore need two versions of both prediction and correction: one
    # doing *only* predict/correct and one doing *only* plot/adjust. The
    # original scheme is correct/adjust -> predict/plot. Modified to allow
    # periodic boundaries: correct -> plot -> adjust -> predict.
    #
    # One could in principle add two extra actions, mappings, adapters,
    # grid functions, etc, but an (ugly) hack minimizing the amount of edits
    # of ExaHyPE is to add two booleans to each Solver: PredictCorrect and
    # PlotAdjust. Only do the prediction/correction step if the relevant boolean
    # is True, and only to the plotting/adjusting step when PlotAdjust is True.

    # We put these booleans in Solver.h so that each solver has them.
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/solvers/Solver.h'

    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    # Check if not already present
    matches = [match for match in lines if "PredictCorrect" in match]
    if (len(matches) == 0):
        print('Modifying Solver.h')
        # By default, they are true, so both steps are fused
        add_class_member(lines, 'Solver', 'public',
                        ['  // Flags to split the prediction/correction step\n',
                         '  bool PredictCorrect = true;\n',
                         '  bool PlotAdjust = true;\n'])

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # We need two extra iterations in the Runner:
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/runners/Runner.cpp'

    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    # NOTE: only works with unfused time steps!!!
    body = remove_function_body(lines,
                                '::runOneTimeStepWithThreeSeparateAlgorithmicSteps')

    # Check if not already present
    matches = [match for match in body if "PredictCorrect" in match]
    if (len(matches) == 0):
        print('Modifying Runner.cpp')

        # Correction step without adjustPointSolution
        for i in range(0, len(body)):
            if (body[i].find('repository.switchToUpdateAndReduce()') != -1):
                body[i+2:i+2] = ['  for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {\n', '    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];\n', '    solver->PredictCorrect = true;\n', '    solver->PlotAdjust = true;\n', '  }\n']
                body[i:i] = ['  for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {\n', '    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];\n', '    solver->PredictCorrect = true;\n', '    solver->PlotAdjust = false;\n', '  }\n']
                break;

        for i in range(0, len(body)):
            if (body[i].find('repository.switchToPrediction()') != -1):
                body[i+2:i+2] = ['  for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {\n', '    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];\n', '    solver->PredictCorrect = true;\n', '    solver->PlotAdjust = true;\n', '  }\n']
                body[i:i] = ['  // Periodic boundaries: first plot\n',
                             '  for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {\n'
                             '    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];\n'
                             '    solver->PredictCorrect = false;\n',
                             '    solver->PlotAdjust = true;\n',
                             '  }\n',
                             '  repository.switchToPrediction();\n',
                             '  repository.iterate( exahype::solvers::Solver::PredictionSweeps, communicatePeanoVertices );\n',
                             '  \n',
                             '  // Periodic boundaries: now adjust solution\n',
                             '  repository.switchToUpdateAndReduce();\n',
                             '  repository.iterate( 1, communicatePeanoVertices );\n',
                             '\n',
                             '  for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {\n'
                             '    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];\n'
                             '    solver->PredictCorrect = true;\n',
                             '    solver->PlotAdjust = false;\n',
                             '  }\n']
                break;

    add_function_body(lines, '::runOneTimeStepWithThreeSeparateAlgorithmicSteps', body)

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # Now update the mappings: first do prediction
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/mappings/Prediction.cpp'

    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    # Replace enterCell function
    body = remove_function_body(lines, '::enterCell')

    # Check if not already present
    matches = [match for match in body if "PredictCorrect" in match]
    if (len(matches) == 0):
        print('Modifying Prediction.cpp')

        for i in range(0, len(body)):
            if (body[i].find('for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++)') != -1):
                body[i+4:i+4] = ['    }\n']
                body[i+3] = '  ' + body[i+3]
                body[i+2] = '  ' + body[i+2]
                body[i+1] = '  ' + body[i+1]
                body[i+0] = '  ' + body[i+0]
                body[i:i] = ['    auto* solver = exahype::solvers::RegisteredSolvers[0];\n', '    if (solver->PlotAdjust == true) {\n']
                break;

        for i in range(0, len(body)):
            if (body[i].find('exahype::mappings::Prediction::performPredictionOrProlongate(') != -1):
                body[i+5:i+5] = ['  }\n']
                body[i+4] = '  ' + body[i+4]
                body[i+3] = '  ' + body[i+3]
                body[i+2] = '  ' + body[i+2]
                body[i+1] = '  ' + body[i+1]
                body[i+0] = '  ' + body[i+0]
                body[i:i] = ['  auto* solver = exahype::solvers::RegisteredSolvers[0];\n', '  if (solver->PredictCorrect == true) {\n']
                break;

    add_function_body(lines, '::enterCell', body)

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    # Now update correction step, which is done in ADERDG
    fname = repo_dir + 'ExaHyPE-Engine/ExaHyPE/exahype/solvers/ADERDGSolver.cpp'

    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    # Replace correction function
    body = remove_function_body(lines, '::correction')

    # Check if not already present
    matches = [match for match in body if "PredictCorrect" in match]
    if (len(matches) == 0):
        print('Modifying ADERDGSolver.cpp')
        for i in range(0, len(body)):
            if (body[i].find('surfaceIntegral(cellDescription,boundaryMarkers,addSurfaceIntegralResultToUpdate)') != -1):
                body[i+4] = '  ' + body[i+4]
                body[i+4:i+4] = ['  }\n', '  if (PlotAdjust == true)\n']
                body[i+3] = '  ' + body[i+3]
                body[i+2] = '  ' + body[i+2]
                body[i+1] = '  ' + body[i+1]
                body[i+0] = '  ' + body[i+0]
                body[i:i] = ['  if (PredictCorrect == true) {\n']
                break;

    add_function_body(lines, '::correction', body)

    f = open(fname, "w")
    f.writelines(lines)
    f.close()
