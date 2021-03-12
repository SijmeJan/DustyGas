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


parser = argparse.ArgumentParser(description='Write flux and eigenvalue cpp')

parser.add_argument('infile',
                    help='ExaHyPE file')
args = parser.parse_args()

f = open(args.infile, "r")
lines = f.readlines()
f.close()

n_vars = 0
output_dir = None
project_name = None
solver_name = None

for line in lines:
    if (line.find('exahype-project ') != -1):
        project_name = line.lstrip().split()[-1]
    if (line.find('solver ') != -1):
        solver_name = line.lstrip().split()[-1]
    if (line.find('variables const') != -1):
        n_vars = int(line.lstrip().split()[-1])
    if (line.find('output-directory') != -1):
        output_dir = line.lstrip().split()[-1]

if (n_vars % 4 != 0):
    print("Error: number of variables has to be divisible by 4!")
    exit(1)

n_dust = int(n_vars/4 - 1)

output_dir = os.path.dirname(os.path.abspath(args.infile)) \
  + '/' + output_dir + '/'

source_file = output_dir + solver_name + '.cpp'

f = open(source_file, "r")
lines = f.readlines()
f.close()

write_eigenvalues(lines, n_dust, 1.0)
write_flux(lines, n_dust, 1.0)

f = open(source_file, "w")
f.writelines(lines)
f.close()
