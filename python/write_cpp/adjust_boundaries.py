import argparse
import os

from scipy.optimize import fsolve

from exahype_parameters import ExaHyPE_parameters
from common import replace_with_indent

def guess_mesh(domain_size, cell_size):
    dx = domain_size
    nx = 1
    while (dx > cell_size):
        nx = 3*nx + 4
        dx = domain_size/nx

    return nx

def periodic_domain_size(domain_size, n_ghost, cell_size):
    # We want a domain that is periodic with domain_size.
    # New size found from: size = domain_size + 2*n_ghost*dx
    f = lambda x: (x - 2*n_ghost*x/guess_mesh(x, cell_size) -
                   domain_size[0])
    size_x_want = fsolve(f, domain_size[0])[0]
    f = lambda x: (x - 2*n_ghost*x/guess_mesh(x, cell_size) -
                   domain_size[1])
    size_y_want = fsolve(f, domain_size[1])[0]

    # Determine mesh number of cells in x and y
    nx = guess_mesh(size_x_want, cell_size)
    ny = guess_mesh(size_y_want, cell_size)

    # Resolution in x and y
    dx = size_x_want/nx
    dy = size_y_want/ny

    print('Old domain size: {}x{}'.format(domain_size[0], domain_size[1]))
    print('New domain size: {}x{}'.format(size_x_want, size_y_want))
    print('Estimated cell size is {}x{}'.format(dx, dy))

    return [size_x_want, size_y_want]

# Needs single argument; exahype file
parser = argparse.ArgumentParser(description='Modify domain size to allow for ghost cells')
parser.add_argument('infile', help='ExaHyPE file')
args = parser.parse_args()

# Read exahype file
project = ExaHyPE_parameters(args.infile)

# Full path to output directory
output_dir = os.path.dirname(os.path.abspath(args.infile)) \
  + '/' + project.output_dir + '/'

# Number of ghost cells required
n_ghost = 1
if project.solver_type == 'Limiting-ADER-DG':
    n_ghost = 2

print('Adjusting domain size to allow for {} ghost layers'.format(n_ghost))

domain_size = \
  periodic_domain_size([project.size_x, project.size_y],
                       n_ghost, project.cell_size)

# We want a domain that is periodic with project.size_x.
# New size found from: size_x = project.size_x + 2*n_ghost*dx
#f = lambda x: (x - 2*n_ghost*x/guess_mesh(x, project.cell_size) -
#               project.size_x)
#size_x_want = fsolve(f, project.size_x)[0]
#f = lambda x: (x - 2*n_ghost*x/guess_mesh(x, project.cell_size) -
#               project.size_y)
#size_y_want = fsolve(f, project.size_y)[0]

# Determine mesh number of cells in x and y
#nx = guess_mesh(size_x_want, project.cell_size)
#ny = guess_mesh(size_y_want, project.cell_size)

# Resolution in x and y
#dx = size_x_want/nx
#dy = size_y_want/ny

#print('Old domain size: {}x{}'.format(project.size_x, project.size_y))
#print('New domain size: {}x{}'.format(size_x_want, size_y_want))
#print('Estimated cell size is {}x{}'.format(dx, dy))

# Adjust exahype file for new domain size
f = open(args.infile, "r")
lines = f.readlines()
f.close()

repl = 'width = {}, {}\n'.format(domain_size[0], domain_size[1])
for i in range(0, len(lines)):
    lines[i] = replace_with_indent(lines[i], 'width = ', repl)

f = open(args.infile, "w")
f.writelines(lines)
f.close()
