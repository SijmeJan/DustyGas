import argparse
import os

from common import replace_with_indent
from eigenvalues import write_eigenvalues
from flux import write_flux
from source import write_source
from initial import write_initial
from boundary import write_boundary, write_boundary_h, write_solver_h
from boundary import write_solver_set_periodic

# Needs single argument; exahype file
parser = argparse.ArgumentParser(description='Write flux and eigenvalue cpp')
parser.add_argument('infile',
                    help='ExaHyPE file')
args = parser.parse_args()

# Read exahype file
f = open(args.infile, "r")
lines = f.readlines()
f.close()

# Parameters to be read from exahype file
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

# Unfortunately, vanilla ExaHyPE does not allow for periodic boundaries.
# A hack that does not require modifying the ExaHyPE core is to make use
# of a plotter that does not generate any output.
use_periodic_boundaries = False

# Find plotter dealing with periodic boundaries
for i in range(0, len(lines)):
    if (lines[i].find('variables const = 0') != -1):
        boundary_name = lines[i-1].lstrip().split()[-1]
        use_periodic_boundaries = True

# Isothermal gas = 4 equations; each dust component adds 4
if (n_vars % 4 != 0):
    print("Error: number of variables has to be divisible by 4!")
    exit(1)

# Parameters
c = 1.0           # gas sound speed
q = 1.5           # non-dimensional shear rate
Stokes = [0.1]    # list of Stokes numbers

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
write_eigenvalues(lines, n_dust, c)
write_flux(lines, n_dust, c)
write_source(lines, n_dust, q, Stokes)
write_initial(lines, n_dust, q, Stokes)

# Write to file
f = open(source_file, "w")
f.writelines(lines)
f.close()

################################
# Implement periodic boundaries
################################
if (use_periodic_boundaries == True):
    print("Implementing periodic boundaries...")

    # Edit solver header file to declare boundary array
    source_file = output_dir + solver_name + '.h'

    f = open(source_file, "r")
    lines = f.readlines()
    f.close()

    write_solver_h(lines)

    f = open(source_file, "w")
    f.writelines(lines)
    f.close()

    # Edit boundary cpp to fill boundary array
    source_file = output_dir + boundary_name + '.cpp'

    f = open(source_file, "r")
    lines = f.readlines()
    f.close()

    write_boundary(lines, n_vars, patch_size,
                   [offset_x, offset_y], [size_x, size_y],
                   solver_name)

    f = open(source_file, "w")
    f.writelines(lines)
    f.close()

    # Edit boundary header to declare pointer to boundary array
    source_file = output_dir + boundary_name + '.h'

    f = open(source_file, "r")
    lines = f.readlines()
    f.close()

    write_boundary_h(lines)

    f = open(source_file, "w")
    f.writelines(lines)
    f.close()


    # Edit solver cpp file to set periodic boundaries
    source_file = output_dir + solver_name + '.cpp'

    f = open(source_file, "r")
    lines = f.readlines()
    f.close()

    write_solver_set_periodic(lines, n_vars)

    f = open(source_file, "w")
    f.writelines(lines)
    f.close()
