from common import replace_with_indent, remove_function_body, add_function_body, add_class_member

def write_outflow_boundary(lines, n_vars, solver):
    # Remove current function body
    remove_function_body(lines, 'boundaryValues')

    periodic = ['  // Outflow boundaries, unimportant if periodic\n']

    if (solver == 'ADER-DG'):
        for n in range(0, n_vars):
            periodic.extend(['  stateOut[{}] = stateIn[{}];\n'.format(n, n)])
        for n in range(0, n_vars):
            periodic.extend(['  fluxOut[{}] = fluxIn[{}];\n'.format(n, n)])
    if (solver == 'Finite-Volumes'):
        for n in range(0, n_vars):
            periodic.extend(['  stateOutside[{}] = stateInside[{}];\n'.format(n, n)])

    add_function_body(lines, 'boundaryValues', periodic)

def write_periodic_functions(n_vars, order, offset, size, output_dir, solver_name):
    fname = output_dir + solver_name + '.cpp'

    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    # Return if function already modified
    matches = [match for match in lines if "int n_cell_x = (int)" in match]
    if len(matches) > 0:
        return

    remove_function_body(lines, '::PlotPeriodic')

    body = \
      ['  //std::cout << "PLOT PERIODIC " << pos[0] << " " << pos[1] << std::endl;\n',
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
      '  // Make sure vector is of the correct size, and set elements to zero.\n',
       '  // Note that this should happen only the first time this function is called.\n',
       '  int n_bound_cells = 2*(n_cell_x + n_cell_y);\n',
       '  int n_send_per_cell = {};\n'.format(n_vars*(order+1)*(order+1)),
       '  if (boundaryVector_local.size() != n_send_per_cell*n_bound_cells) {\n',
       '    boundaryVector_local.resize(n_send_per_cell*n_bound_cells);\n',
       '    std::fill(boundaryVector_local.begin(), boundaryVector_local.end(), 0.0);\n',
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
       '      boundaryVector_local[arr_index + n] = Q[n];\n',
       '  }\n'
       '  if (cell_y == n_cell_y - 2) {\n',
       '    indx = n_cell_x + cell_x;\n',
       '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
       '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
       '      boundaryVector_local[arr_index + n] = Q[n];\n',
       '  }\n'
       '  if (cell_x == 1) {\n',
       '    indx = 2*n_cell_x + cell_y;\n',
       '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
       '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
       '      boundaryVector_local[arr_index + n] = Q[n];\n',
       '  }\n'
       '  if (cell_x == n_cell_x - 2) {\n',
       '    indx = 2*n_cell_x + n_cell_y + cell_y;\n',
       '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
       '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
       '      boundaryVector_local[arr_index + n] = Q[n];\n',
       '  }\n']

    add_function_body(lines, '::PlotPeriodic', body)

    remove_function_body(lines, '::SendPeriodic')

    body = \
      ['  //std::cout << "SENDING PERIODIC" << std::endl;\n',
       '#ifdef Parallel\n',
       '  // Should only happen the first time of call\n',
       '  if (boundaryVector.size() == 0) {\n',
       '    if (tarch::parallel::Node::getInstance().getNumberOfNodes() > 1) {\n',
       '      // Main rank does not know the boundary size\n',
       '      int s_loc = boundaryVector_local.size();\n',
       '      int s_max = 0;\n',
       '      MPI_Reduce(&s_loc, &s_max, 1, MPI_INT, MPI_MAX, 0,\n',
       '                 tarch::parallel::Node::getInstance().getCommunicator());\n',
       '      // Set correct size so we can use Allreduce\n',
       '      if (tarch::parallel::Node::getInstance().isGlobalMaster()) {\n',
       '        boundaryVector_local.resize(s_max);\n',
       '        std::fill(boundaryVector_local.begin(), boundaryVector_local.end(), 0.0);\n',
       '        boundaryVector.resize(s_max);\n',
       '      } else {\n',
       '        boundaryVector.resize(boundaryVector_local.size());\n',
       '      }\n',
       '    } else {\n',
       '      boundaryVector.resize(boundaryVector_local.size());\n',
       '    }\n',
       '  }\n',
       '  //std::cout << "Allreduce with rank " << tarch::parallel::Node::getInstance().getRank() << " and communicator " << tarch::parallel::Node::getInstance().getCommunicator() << std::endl;\n',
       '  MPI_Allreduce(&(boundaryVector_local[0]),\n',
       '                &(boundaryVector[0]),\n',
       '                boundaryVector_local.size(),\n',
       '                MPI_DOUBLE,\n',
       '                MPI_SUM,\n',
       '                tarch::parallel::Node::getInstance().getCommunicator());\n',
       '#else\n',
       '  boundaryVector = boundaryVector_local;\n',
       '#endif\n']

    add_function_body(lines, '::SendPeriodic', body)

    remove_function_body(lines, '::AdjustPeriodic')

    body = \
      ['  //std::cout << "ADJUSTPERIODIC " << pos[0] << " " << pos[1] << std::endl;\n',
       '  // Hack: number of cells in x and y\n',
       '  int n_cell_x = (int) round({}/sizeOfPatch[0]);\n'.format(size[0]),
       '  int n_cell_y = (int) round({}/sizeOfPatch[1]);\n'.format(size[1]),
       '\n',
       '  // Number of cells to left and bottom\n',
       '  int cell_x = (int) ((offsetOfPatch[0] + 0.5*sizeOfPatch[0])/sizeOfPatch[0]);\n',
       '  int cell_y = (int) ((offsetOfPatch[1] + 0.5*sizeOfPatch[1])/sizeOfPatch[1]);\n',
       '\n',
       '  int n_send_per_cell = {};\n'.format(n_vars*(order+1)*(order+1)),
       '  int indx = -1;\n',
       '\n',
       '  if (cell_y == n_cell_y - 1) {\n',
       '    indx = cell_x;\n',
       '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
       '    //std::cout << "ADUSTING: " << cell_x << " " << cell_y << " " << indx << " " << arr_index << " " << boundaryVector.size() << std::endl;\n',
       '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
       '      Q[n] = boundaryVector[arr_index + n];\n',
       '  }\n'
       '  if (cell_y == 0) {\n',
       '    indx = n_cell_x + cell_x;\n',
       '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
       '    //std::cout << "ADUSTING: " << cell_x << " " << cell_y << " " << indx << " " << arr_index << " " << boundaryVector.size() << std::endl;\n',
       '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
       '      Q[n] = boundaryVector[arr_index + n];\n',
       '  }\n'
       '  if (cell_x == n_cell_x - 1) {\n',
       '    indx = 2*n_cell_x + cell_y;\n',
       '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
       '    //std::cout << "ADUSTING: " << cell_x << " " << cell_y << " " << indx << " " << arr_index << " " << boundaryVector.size() << std::endl;\n',
       '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
       '      Q[n] = boundaryVector[arr_index + n];\n',
       '  }\n'
       '  if (cell_x == 0) {\n',
       '    indx = 2*n_cell_x + n_cell_y + cell_y;\n',
       '    int arr_index = n_send_per_cell*indx + ({}*pos[1] + pos[0])*{};\n'.format(order + 1, n_vars),
       '    //std::cout << "ADUSTING: " << cell_x << " " << cell_y << " " << indx << " " << arr_index << " " << boundaryVector.size() << std::endl;\n',
       '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
       '      Q[n] = boundaryVector[arr_index + n];\n',
       '  }\n'
      ]

    add_function_body(lines, '::AdjustPeriodic', body)

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

def write_periodic_dummies(base_filename, solver_name):
    fname = base_filename + '.h'

    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    # Return if function already present
    matches = [match for match in lines if "void AdjustPeriodic" in match]
    if len(matches) > 0:
        return

    for i in range(1, len(lines)):
        if (lines[-i].find('};') != -1):
            lines[-i:-i] = \
              ['\n',
               '  void AdjustPeriodic(\n',
               '    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
               '    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
               '    const tarch::la::Vector<DIMENSIONS, int>& pos,\n',
               '    double* Q) override;\n',
               '  void PlotPeriodic(\n',
               '    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
               '    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
               '    const tarch::la::Vector<DIMENSIONS, int>& pos,\n',
               '    double* const Q) override;\n',
               '  void SendPeriodic() override;\n',
               '\n',
               '  std::vector<double> boundaryVector;\n',
               '  std::vector<double> boundaryVector_local;\n',
               '  //double *global_dx;\n',
               '  //int *global_n;\n']
            break;

    # Make sure to include vector
    boundary = ['#include <vector>\n']

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i:i] = boundary
            break;

    f = open(fname, "w")
    f.writelines(lines)
    f.close()

    fname = base_filename + '.cpp'

    f = open(fname, "r")
    lines = f.readlines()
    f.close()

    lines[len(lines):len(lines)] = \
      ['\n\n',
       'void {}::{}::AdjustPeriodic(\n'.format(solver_name[:-6], solver_name),
       '    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
       '    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
       '    const tarch::la::Vector<DIMENSIONS, int>& pos,\n',
       '    double* Q) {\n'
       '  // Nop, since no periodic boundaries requested\n',
       '}\n\n',
       'void {}::{}::PlotPeriodic(\n'.format(solver_name[:-6], solver_name),
       '    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,\n',
       '    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,\n',
       '    const tarch::la::Vector<DIMENSIONS, int>& pos,\n',
       '    double* const Q) {\n'
       '  // Nop, since no periodic boundaries requested\n',
       '}\n\n',
       'void {}::{}::SendPeriodic() {{\n'.format(solver_name[:-6], solver_name),
       '  // Nop, since no periodic boundaries requested\n',
       '}\n\n']

    f = open(fname, "w")
    f.writelines(lines)
    f.close()
