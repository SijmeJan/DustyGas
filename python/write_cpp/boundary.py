from common import replace_with_indent, remove_function_body, add_function_body

def write_boundary(lines, n_vars, patch_size, offset, size, solver_name):
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
                '  // Number of patches to left and bottom\n',
                '  int patch_x = (int) ((offsetOfPatch[0] + 0.5*sizeOfPatch[0])/sizeOfPatch[0]);\n',
                '  int patch_y = (int) ((offsetOfPatch[1] + 0.5*sizeOfPatch[1])/sizeOfPatch[1]);\n',
                '\n',
                '  int indx = -1;\n',
                '\n',
                '  if (x[1] - 0.5*sizeOfPatch[1] < {} && pos[1] == 0) {{\n'.format(y_bound[0]),
                '    // Bottom boundary\n',
                '    indx = patch_x*{} + pos[0];\n'.format(patch_size),
                '  }\n',
                '  if (x[1] + 0.5*sizeOfPatch[1] > {} && pos[1] == {}) {{\n'.format(y_bound[1], patch_size-1),
                '    // Top boundary\n',
                '    indx = (n_patch_x + patch_x)*{} + pos[0];\n'.format(patch_size),
                '  }\n',
                '  if (x[0] - 0.5*sizeOfPatch[0] < {} && pos[0] == 0) {{\n'.format(x_bound[0]),
                '    // Left boundary\n',
                '    indx = (2*n_patch_x + patch_y)*{} + pos[1];\n'.format(patch_size),
                '  }\n',
                '  if (x[0] + 0.5*sizeOfPatch[0] > {} && pos[0] == {}) {{\n'.format(x_bound[1], patch_size-1),
                '    // Right boundary\n',
                '    indx = (2*n_patch_x + n_patch_y + patch_y)*{} + pos[1];\n'.format(patch_size),
                '  }\n',
                '\n',
                '  if (indx > -1) {\n',
                '    int arr_index = {}*indx;\n'.format(n_vars),
                '    // Resize if necessary\n',
                '    if (arr_index >= (*boundaryValues).size())\n',
                '      (*boundaryValues).resize(arr_index + {});\n'.format(n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      (*boundaryValues)[arr_index + n] = Q[n];\n',
                '  }\n']

    add_function_body(lines, 'mapQuantities', boundary)

    # Include solver header
    solver = ['#include "{}.h"\n'.format(solver_name)]

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i+1:i+1] = solver
            break;

    solver = ['  boundaryValues = &(solver.periodicBoundaryValues);\n']

    for i in range(0, len(lines)):
        if (lines[i].find('::') != -1):
            lines[i+1:i+1] = solver
            break;

def write_boundary_h(lines):
    boundary = [' private:\n',
                '  std::vector<double> *boundaryValues;\n']

    # Do nothing if line already present
    for i in range(0, len(lines)):
        if (lines[i].find(boundary[1]) != -1):
            return

    for i in range(0, len(lines)):
        if (lines[i].find('public:') != -1):
            lines[i:i] = boundary
            break;

    boundary = ['#include <vector>\n']

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i:i] = boundary
            break;

def write_solver_h(lines):
    solver = ['    std::vector<double> periodicBoundaryValues;\n']

    # Do nothing if line already present
    for i in range(0, len(lines)):
        if (lines[i].find(solver[0]) != -1):
            return

    for i in range(0, len(lines)):
        if (lines[i].find('};') != -1):
            lines[i:i] = solver
            break;

    solver = ['#include <vector>\n']

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i:i] = solver
            break;

def write_solver_set_periodic(lines):
    # Save current function body
    body = remove_function_body(lines, 'adjustSolution')

    periodic = ['  if (t > 0.0) {\n',
                '    std::cout << x[0] << " " << x[1] << std::endl;\n',
                '  }\n']

    body.extend(periodic)

    add_function_body(lines, 'adjustSolution', body)
