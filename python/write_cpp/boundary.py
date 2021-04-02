from common import replace_with_indent, remove_function_body, add_function_body, add_class_member

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
                '    std::cout << "Mapping quantities at x = " << x[0] << ", y = " << x[1] << ", position " << pos[0] << " " << pos[1] << ", index = " << indx << ", arr_index = " << arr_index << ", cell index " << cell_x << " " << cell_y << std::endl;\n',
                '  }\n'
                '  if (cell_y == global_n[1] - 2) {\n',
                '    indx = n_cell_x + cell_x;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '    std::cout << "Mapping quantities at x = " << x[0] << ", y = " << x[1] << ", position " << pos[0] << " " << pos[1] << ", index = " << indx << ", arr_index = " << arr_index << ", cell index " << cell_x << " " << cell_y << std::endl;\n',
                '  }\n'
                '  if (cell_x == 1) {\n',
                '    indx = 2*n_cell_x + cell_y;\n',
                '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      boundaryValues_local[arr_index + n] = Q[n];\n',
                '    std::cout << "Mapping quantities at x = " << x[0] << ", y = " << x[1] << ", position " << pos[0] << " " << pos[1] << ", index = " << indx << ", arr_index = " << arr_index << ", cell index " << cell_x << " " << cell_y << std::endl;\n',
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

    periodic = ['  if (t > 0.0) {\n',
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
                '      if (i == 0 && j == 0) std::cout << "Adjusting point solution at x = " << x[0] << ", y = " << x[1] << ", cell index " << i << " " << j << ", index = " << indx << ", arr_index = " << arr_index << " at t = " << t << std::endl;\n',
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
    remove_function_body(lines, 'boundaryValues')

    periodic = ['  // Outflow boundaries, unimportant since periodic anyway\n']

    for n in range(0, n_vars):
        periodic.extend(['  stateOut[{}] = stateIn[{}];\n'.format(n, n)])
    for n in range(0, n_vars):
        periodic.extend(['  fluxOut[{}] = fluxIn[{}];\n'.format(n, n)])

    add_function_body(lines, 'boundaryValues', periodic)

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


def correction_boundary_hack(repo_dir):
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
                body[i:i] = ['  for (int solverNumber=0; solverNumber<static_cast<int>(solvers::RegisteredSolvers.size()); solverNumber++) {\n', '    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];\n', '    solver->PredictCorrect = true;\n', '    solver->PlotAdjust = false;\n', '  }\n']
                break;

        for i in range(0, len(body)):
            if (body[i].find('repository.switchToPrediction()') != -1):
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
                body[i:i] = ['  if (solver->PredictCorrect == true) {\n']
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
