import argparse
import os

def replace_with_indent(line, search, replace):
    ret = line

    if (line.find(search) != -1):
        n_white  = len(line) - len(line.lstrip())
        ret = replace.rjust(len(replace) + n_white)

    return ret

def write_eigenvalues(lines, n_dust, c):
    for i in range(0, len(lines)):
        lines[i] = replace_with_indent(lines[i],
                                      'lambda[0] = ',
                                      'lambda[0] = Q[dIndex + 1]/Q[0] - {};\n'.format(c))
        lines[i] = replace_with_indent(lines[i],
                                      'lambda[1] = ',
                                      'lambda[1] = Q[dIndex + 1]/Q[0];\n')
        lines[i] = replace_with_indent(lines[i],
                                      'lambda[2] = ',
                                      'lambda[2] = Q[dIndex + 1]/Q[0];\n')
        lines[i] = replace_with_indent(lines[i],
                                      'lambda[3] = ',
                                      'lambda[3] = Q[dIndex + 1]/Q[0] + {};\n'.format(c))
        for n in range(0, n_dust):
            for j in range(0, 4):
                lines[i] = replace_with_indent(lines[i],
                                              'lambda[{}] = '.format(4*n + 4 + j),
                                              'lambda[{}] = Q[dIndex + {}]/Q[{}];\n'.format(4*n + 4 + j, 5 + 4*n, 4 + 4*n))

def write_flux(lines, n_dust, c):
    for i in range(0, len(lines)):
        # First direction
        lines[i] = replace_with_indent(lines[i],
                                      'F[0][0] = ',
                                      'F[0][0] = Q[1];\n')
        lines[i] = replace_with_indent(lines[i],
                                      'F[0][1] = ',
                                      'F[0][1] = Q[1]*Q[1]/Q[0] + {}*Q[0];\n'.format(c*c))
        lines[i] = replace_with_indent(lines[i],
                                      'F[0][2] = ',
                                      'F[0][2] = Q[1]*Q[2]/Q[0];\n')
        lines[i] = replace_with_indent(lines[i],
                                      'F[0][3] = ',
                                      'F[0][3] = Q[1]*Q[3]/Q[0];\n')

        # Second direction
        lines[i] = replace_with_indent(lines[i],
                                      'F[1][0] = ',
                                      'F[1][0] = Q[2];\n')
        lines[i] = replace_with_indent(lines[i],
                                      'F[1][1] = ',
                                      'F[1][1] = Q[2]*Q[1]/Q[0];\n')
        lines[i] = replace_with_indent(lines[i],
                                      'F[1][2] = ',
                                      'F[1][2] = Q[2]*Q[2]/Q[0] + {}*Q[0];\n'.format(c*c))
        lines[i] = replace_with_indent(lines[i],
                                      'F[1][3] = ',
                                      'F[1][3] = Q[2]*Q[3]/Q[0];\n')


        for n in range(0, n_dust):
            lines[i] = replace_with_indent(lines[i],
                                           'F[0][{}] = '.format(4*n + 4),
                                           'F[0][{}] = Q[{}];\n'.format(4*n + 4, 5 + 4*n))
            lines[i] = replace_with_indent(lines[i],
                                           'F[0][{}] = '.format(4*n + 5),
                                           'F[0][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 5, 4*n + 5, 4*n + 5, 4*n + 4))
            lines[i] = replace_with_indent(lines[i],
                                           'F[0][{}] = '.format(4*n + 6),
                                           'F[0][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 6, 4*n + 5, 4*n + 6, 4*n + 4))
            lines[i] = replace_with_indent(lines[i],
                                           'F[0][{}] = '.format(4*n + 7),
                                           'F[0][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 7, 4*n + 5, 4*n + 7, 4*n + 4))



            lines[i] = replace_with_indent(lines[i],
                                           'F[1][{}] = '.format(4*n + 4),
                                           'F[1][{}] = Q[{}];\n'.format(4*n + 4, 4*n + 6))
            lines[i] = replace_with_indent(lines[i],
                                           'F[1][{}] = '.format(4*n + 5),
                                           'F[1][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 5, 4*n + 6, 4*n + 5, 4*n + 4))
            lines[i] = replace_with_indent(lines[i],
                                           'F[1][{}] = '.format(4*n + 6),
                                           'F[1][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 6, 4*n + 6, 4*n + 6, 4*n + 4))
            lines[i] = replace_with_indent(lines[i],
                                           'F[1][{}] = '.format(4*n + 7),
                                           'F[1][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 7, 4*n + 6, 4*n + 7, 4*n + 4))

def write_source(lines, n_dust, q, Stokes):

    for i in range(0, len(lines)):
        # Source x: eta + Coriolis
        lines[i] = replace_with_indent(lines[i],
                                      'S[1] = ',
                                      'S[1] = 2*Q[0] + 2*Q[3];\n')
        # Source y: Coriolis
        lines[i] = replace_with_indent(lines[i],
                                      'S[3] = ',
                                      'S[3] = ({} - 2)*Q[1];\n'.format(q))

        for n in range(0, n_dust):
            # Source x: Coriolis + drag
            lines[i] = replace_with_indent(lines[i],
                                           'S[{}] = '.format(4*n + 5),
                                           'S[{}] = 2*Q[{}] - (Q[{}] - Q[{}]*Q[1]/Q[0])/{};\n'.format(4*n + 5, 4*n + 7, 4*n + 5, 4*n + 4, Stokes[n]))

            # Source z: drag
            lines[i] = replace_with_indent(lines[i],
                                           'S[{}] = '.format(4*n + 6),
                                           'S[{}] = - (Q[{}] - Q[{}]*Q[1]/Q[0])/{};\n'.format(4*n + 6, 4*n + 6, 4*n + 4, Stokes[n]))

            # Source y: Coriolis + drag
            lines[i] = replace_with_indent(lines[i],
                                          'S[{}] = '.format(4*n + 7),
                                          'S[{}] = ({} - 2)*Q[{}] - (Q[{}] - Q[{}]*Q[1]/Q[0])/{};\n'.format(4*n + 7, q, 4*n + 5, 4*n + 7, 4*n + 4, Stokes[n]))

    gas_x_drag = []
    gas_y_drag = []
    gas_z_drag = []
    for n in range(0, n_dust):
        gas_x_drag.append('  S[1] += (Q[{}] -Q[{}]*Q[1]/Q[0])/{};\n'.format(4*n + 5, 4*n + 4, Stokes[n]))
        gas_z_drag.append('  S[2] += (Q[{}] -Q[{}]*Q[1]/Q[0])/{};\n'.format(4*n + 6, 4*n + 4, Stokes[n]))
        gas_y_drag.append('  S[3] += (Q[{}] -Q[{}]*Q[1]/Q[0])/{};\n'.format(4*n + 7, 4*n + 4, Stokes[n]))

    for i in range(0, len(lines)):
        if (lines[i].find('S[1] = ') != -1):
            lines[i+1:i+1] = gas_x_drag
    for i in range(0, len(lines)):
        if (lines[i].find('S[2] = ') != -1):
            lines[i+1:i+1] = gas_z_drag
    for i in range(0, len(lines)):
        if (lines[i].find('S[3] = ') != -1):
            lines[i+1:i+1] = gas_y_drag

def write_initial(lines, n_dust, q, Stokes):
    for i in range(0, len(lines)):
        for n in range(0, 4*n_dust + 4):
            lines[i] = replace_with_indent(lines[i],
                                           'Q[{}] = 0.0'.format(n),
                                           '//Q[{}] = 0.0;\n'.format(n))

    initial = ['  if (t == 0.0) {\n',
               '    Q[0] = 1.0;\n',
               '    Q[1] = 0.0;\n',
               '    Q[2] = 0.0;\n',
               '    Q[3] = 0.0;\n',
               '    if (sqrt((x[0] - 0.7)*(x[0] - 0.7) + (x[1] - 0.7)*(x[1] - 0.7)) < 0.2) Q[0] = 2.0;\n',
               '  }\n']

    for i in range(0, len(lines)):
        if (lines[i].find('//Q[0] = ') != -1):
            lines[i:i] = initial
            break;

def write_boundary(lines, n_vars, patch_size, offset, size):
    for i in range(0, len(lines)):
        lines[i] = replace_with_indent(lines[i],
                                      'const int writtenUnknowns = 0;',
                                      'assertion(outputQuantities==nullptr);\n')

    for i in range(0, len(lines)):
        if (lines[i].find('assertion') != -1):
            lines.pop(i+1)
            lines.pop(i+1)
            lines.pop(i+1)
            break;

    x_bound = [offset[0], offset[0] + size[0]]
    y_bound = [offset[1], offset[1] + size[1]]



    boundary = ['  // Fill a boundary array for setting periodic boundaries in 2D, non-AMR runs.\n',
                '  // If FV mesh = nx times ny, the first nx entries correspond to the bottom boundary.\n',
                '  // The second nx entries correspond to the top boundary.\n',
                '  // The next ny entries correspond to the left boundary.\n',
                '  // The next ny entries correspont to the right boundary.\n',
                '\n',
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
                '    if (arr_index >= solver.periodicBoundaryValues.size())\n',
                '      solver.periodicBoundaryValues.resize(arr_index + {});\n'.format(n_vars),
                '    for (int n = 0; n < {}; n++)\n'.format(n_vars),
                '      solver.periodicBoundaryValues[arr_index + n] = Q[n];\n',
                '  }\n']

    for i in range(0, len(lines)):
        if (lines[i].find('assertion') != -1):
            lines[i+1:i+1] = boundary
            break;

def write_boundary_h(lines):
    boundary = [' private:\n',
                '  std::vector<double> boundaryValues;\n']

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
    solver = ['  std::vector<double> periodicBoundaryValues;\n']

    for i in range(0, len(lines)):
        if (lines[i].find('};') != -1):
            lines[i:i] = solver
            break;

    solver = ['#include <vector>\n']

    for i in range(0, len(lines)):
        if (lines[i].find('#include') != -1):
            lines[i:i] = solver
            break;

parser = argparse.ArgumentParser(description='Write flux and eigenvalue cpp')

parser.add_argument('infile',
                    help='ExaHyPE file')
args = parser.parse_args()

f = open(args.infile, "r")
lines = f.readlines()
f.close()

n_vars = 0
patch_size = 0
output_dir = None
project_name = None
solver_name = None
boundary_name = None
found_nvar = False
offset_x = None
offset_y = None
size_x = None
size_y = None


for line in lines:
    if (line.find('exahype-project ') != -1):
        project_name = line.lstrip().split()[-1]
    if (line.find('solver ') != -1):
        solver_name = line.lstrip().split()[-1]
    if (line.find('variables const') != -1 and found_nvar == False):
        n_vars = int(line.lstrip().split()[-1])
        found_nvar = True
    if (line.find('output-directory') != -1):
        output_dir = line.lstrip().split()[-1]
    if (line.find('width') != -1):
        size_x = float(line.lstrip().split()[-2][:-1])
        size_y = float(line.lstrip().split()[-1])
    if (line.find('offset') != -1):
        offset_x = float(line.lstrip().split()[-2][:-1])
        offset_y = float(line.lstrip().split()[-1])
    if (line.find('patch-size const') != -1):
        patch_size = int(line.lstrip().split()[-1])

for i in range(0, len(lines)):
    if (lines[i].find('variables const = 0') != -1):
        boundary_name = lines[i-1].lstrip().split()[-1]

if (n_vars % 4 != 0):
    print("Error: number of variables has to be divisible by 4!")
    exit(1)

c = 1.0
q = 1.5
Stokes = [0.1]

n_dust = int(n_vars/4 - 1)

output_dir = os.path.dirname(os.path.abspath(args.infile)) \
  + '/' + output_dir + '/'

source_file = output_dir + solver_name + '.cpp'

f = open(source_file, "r")
lines = f.readlines()
f.close()

write_eigenvalues(lines, n_dust, c)
write_flux(lines, n_dust, c)
write_source(lines, n_dust, q, Stokes)
write_initial(lines, n_dust, q, Stokes)

f = open(source_file, "w")
f.writelines(lines)
f.close()


# PERIODIC BOUNDARIES HACK
#
# Naive (but robust?) implementation:
# - global array of boundary values for *whole* grid
# - all workers send boundary data to array (MPI_Allgather "plotting")
# - adjust solution from array
# Advantage:
# - probably works with loadbalancing
# - no need to touch ExaHyPE core
# Disadvantage:
# - lot of communication
# - no AMR (this would in any case require a serious hack)
#
# Optimal implementation:
# - set rank of boundary vertices and let Peano handle the rest...
# Advantage: minimal communication
# Disadvantage: hack into Peano, no AMR, no loadbalancing (?)

source_file = output_dir + solver_name + '.h'

f = open(source_file, "r")
lines = f.readlines()
f.close()

write_solver_h(lines)

f = open(source_file, "w")
f.writelines(lines)
f.close()


source_file = output_dir + boundary_name + '.cpp'

f = open(source_file, "r")
lines = f.readlines()
f.close()

write_boundary(lines, n_vars, patch_size,
               [offset_x, offset_y], [size_x, size_y])

f = open(source_file, "w")
f.writelines(lines)
f.close()

#source_file = output_dir + boundary_name + '.h'

#f = open(source_file, "r")
#lines = f.readlines()
#f.close()

#write_boundary_h(lines)

#f = open(source_file, "w")
#f.writelines(lines)
#f.close()
