import argparse
import os
from scipy.optimize import fsolve

from common import replace_with_indent
from eigenvalues import write_eigenvalues
from flux import write_flux
from source import write_source
from initial import write_initial
from boundary import write_boundary, write_boundary_h, write_solver_h
from boundary import write_solver_set_periodic
from boundary import write_outflow_boundary
from boundary import write_periodic_functions
from boundary import write_periodic_dummies
from plot import write_plotter_gas_velocity
from plot import write_plotter_dust_velocity
from plot import write_plotter_dust_density
from periodic_exahype import allow_periodic

def guess_mesh(domain_size, cell_size):
    dx = domain_size
    nx = 1
    while (dx > cell_size):
        nx = 3*nx + 4
        dx = domain_size/nx

    return nx

# Needs single argument; exahype file
parser = argparse.ArgumentParser(description='Write flux and eigenvalue cpp')
parser.add_argument('--periodic', action="store_true", help="use periodic domain")
parser.add_argument('infile',
                    help='ExaHyPE file')
args = parser.parse_args()


# Read exahype file
f = open(args.infile, "r")
lines = f.readlines()
f.close()

# Remove comments between /* and */
start_comment = []
end_comment = []
for i in range(0, len(lines)):
    if (lines[i].find('/*') != -1):
        start_comment.append(i)
    if (lines[i].find('*/') != -1):
        end_comment.append(i)
for n in range(len(start_comment)-1, -1, -1):
    for i in range(0, end_comment[n] - start_comment[n] + 1):
        lines.pop(start_comment[n])

# Parameters to be read from exahype file
n_vars = 0
patch_size = 0
order = None
output_dir = None
project_name = None
solver_name = None
solver_type = None
boundary_name = None
found_nvar = False
offset_x = None
offset_y = None
size_x = None
size_y = None
cell_size = None

for line in lines:
    if (line.find('exahype-project ') != -1):
        project_name = line.lstrip().split()[-1]
    if (line.find('solver ') != -1):
        solver_name = line.lstrip().split()[-1]
        solver_type = line.lstrip().split()[-2]
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
    if (line.find('order const') != -1):
        order = int(line.lstrip().split()[-1])
    if (line.find('maximum-mesh-size') != -1):
        cell_size = float(line.lstrip().split()[-1])

# Main repository directory
repo_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../'

# User signals they want periodic domain
# Check if ExaHyPE is set up to do this
if args.periodic:
    fname = repo_dir + \
      'ExaHyPE-Engine/ExaHyPE/exahype/mappings/PlotPeriodic.cpp'
    if os.path.isfile(fname) != True:
        print('ExaHyPE is not set up for periodic domains.')
        print('Run allow_periodic.sh before creating the project.')
        exit(1)

# Unfortunately, vanilla ExaHyPE does not allow for periodic boundaries.
# Check if periodic boundaries have been enabled?
use_periodic_boundaries = True

# Isothermal gas = 4 equations; each dust component adds 4
if (n_vars % 4 != 0):
    print("Error: number of variables has to be divisible by 4!")
    exit(1)

# Parameters
c = 1.0           # gas sound speed
eta = 0.0         # eta/c

q = 1.5           # non-dimensional shear rate
Stokes = [0.1]    # list of Stokes numbers
mu = 3.0          # dust/gas ratio

# Number of dust components
n_dust = int(n_vars/4 - 1)

# Full path to cpp files
output_dir = os.path.dirname(os.path.abspath(args.infile)) \
  + '/' + output_dir + '/'

#####################################
# Start by modifying the solver file
#####################################
source_file = output_dir + solver_name + '.cpp'

f = open(source_file, "r")
lines = f.readlines()
f.close()

# Modify eigenvalues, fluxes, sources and initial conditions
write_eigenvalues(lines, n_dust, c, solver_type)
write_flux(lines, n_dust, c)
write_source(lines, n_dust, q, Stokes, eta)
write_initial(lines, n_dust, mu, Stokes, eta, solver_type)
write_outflow_boundary(lines, n_vars)

# Write to file
f = open(source_file, "w")
f.writelines(lines)
f.close()

# Hack into Exahype to allow two extra mappings
#allow_periodic(repo_dir)

# Write empty functions (no periodic boundaries by default)
# Can only do this if periodic boundaries have been enabled
write_periodic_dummies(output_dir, solver_name)

################################
# Implement periodic boundaries
################################
if (use_periodic_boundaries == True):
    print("Implementing periodic boundaries...")

    if (solver_type == 'ADER-DG'):
        nx = guess_mesh(size_x, cell_size)
        ny = guess_mesh(size_y, cell_size)

        dx = size_x/nx
        dy = size_y/ny

        f = lambda x: x - 2*x/guess_mesh(x, cell_size) - size_x
        size_x_want = fsolve(f, size_x)[0]
        f = lambda x: x - 2*x/guess_mesh(x, cell_size) - size_y
        size_y_want = fsolve(f, size_y)[0]

        print('Estimated mesh size: {}x{}'.format(nx, ny))
        print('NOTE: periodic domain size will be {}x{}'.format(size_x-2*dx, size_y-2*dy))
        print('If a periodic domain of size {}x{} is needed, change domain size in .exahype file to {}x{}'.format(size_x, size_y, size_x_want, size_y_want))

    # Edit solver header file to declare boundary array
    #source_file = output_dir + solver_name + '.h'

    #f = open(source_file, "r")
    #lines = f.readlines()
    #f.close()

    #write_solver_h(lines)

    #f = open(source_file, "w")
    #f.writelines(lines)
    #f.close()

    # Edit boundary cpp to fill boundary array
    #source_file = output_dir + boundary_name + '.cpp'

    #f = open(source_file, "r")
    #lines = f.readlines()
    #f.close()

    #write_boundary(lines, n_vars, order,
    #               [offset_x, offset_y], [size_x, size_y],
    #               solver_name, boundary_name)

    #f = open(source_file, "w")
    #f.writelines(lines)
    #f.close()

    # Edit boundary header to declare pointer to boundary array
    #source_file = output_dir + boundary_name + '.h'

    #f = open(source_file, "r")
    #lines = f.readlines()
    #f.close()

    #write_boundary_h(lines)

    #f = open(source_file, "w")
    #f.writelines(lines)
    #f.close()


    # Edit solver cpp file to set periodic boundaries
    #source_file = output_dir + solver_name + '.cpp'

    #f = open(source_file, "r")
    #lines = f.readlines()
    #f.close()

    #write_solver_set_periodic(lines, n_vars, order)

    #f = open(source_file, "w")
    #f.writelines(lines)
    #f.close()


    # Add PlotAdjust function to abstract solver class
    write_periodic_functions(n_vars, order,
                             [offset_x, offset_y],
                             [size_x, size_y],
                             output_dir, solver_name)
