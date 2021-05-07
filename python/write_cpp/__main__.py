import argparse
import os
import numpy as np
import configparser

from scipy.optimize import fsolve
from scipy.special import roots_legendre

from eigenvalues import write_eigenvalues
from flux import write_flux
from source import write_source
from initial import write_initial
from boundary import write_outflow_boundary
from boundary import write_periodic_functions
from boundary import write_periodic_dummies
from physical import write_physical
from size_density import SizeDensity
from exahype_parameters import ExaHyPE_parameters

# Needs single argument; exahype file
parser = argparse.ArgumentParser(description='Write flux and eigenvalue cpp')
parser.add_argument('--periodic', action="store_true", help="use periodic domain")
parser.add_argument('infile',
                    help='ExaHyPE file')
args = parser.parse_args()

# Read exahype file
project = ExaHyPE_parameters(args.infile)

# Main repository directory
repo_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../'

# Check if ExaHyPE is set up to do periodic domains
allow_periodic = False
fname = repo_dir + \
  'ExaHyPE-Engine/ExaHyPE/exahype/mappings/PlotPeriodic.cpp'
if os.path.isfile(fname) == True:
    allow_periodic = True

# If we want a periodic domain, make sure it is set up
if args.periodic:
    if allow_periodic == False:
        print('ExaHyPE is not set up for periodic domains.')
        print('Run allow_periodic.sh before creating the project.')
        exit(1)

# Isothermal gas = 4 equations; each dust component adds 4
if (project.n_vars % 4 != 0):
    print("Error: number of variables has to be divisible by 4!")
    exit(1)

# Read dust/gas parameters from input file
ini_file = args.infile.split('.exahype')[0] + '.ini'
# Parse INI file
config = configparser.ConfigParser()
config.read(ini_file)

Kx = float(config['DOMAIN']['Kx'])
Kz = float(config['DOMAIN']['Kz'])
c = float(config['GAS']['c'])
eta = 1.0
q = float(config['GAS']['q'])
Stokes = [float(config['DUST']['Stokes'])]
mu = float(config['DUST']['mu'])

# Number of dust components
n_dust = int(project.n_vars/4 - 1)

# Full path to solver cpp files to write
output_dir = os.path.dirname(os.path.abspath(args.infile)) \
  + '/' + project.output_dir + '/'

# If we are doing a size continuum...
weights = [1.0]
Stokes_range = [0.0999, 0.1]
sigma=None
if n_dust > 1:
    # Get nodes and weights for Gauss-Legendre
    x, w = roots_legendre(n_dust)
    # Adjust for integrals from 0 to 1
    lk = 0.5*(x + 1)
    wk = 0.5*w

    # Work with log(St)
    s_range = np.log(np.asarray(Stokes_range))

    Stokes = np.exp(s_range[0] + (s_range[1] - s_range[0])*lk)
    weights = (s_range[1] - s_range[0])*wk

    # MRN size density
    sigma = SizeDensity(0.5, Stokes_range)

#####################################
# Start by modifying the solver file
#####################################
solver_types = [project.solver_type]
solver_extension = ['']

# For a limiting scheme, need to modify both solvers
if (project.solver_type == 'Limiting-ADER-DG'):
    solver_types = ['Finite-Volumes', 'ADER-DG']
    solver_extension = ['_FV', '_ADERDG']

for s_ext, s_type in zip(solver_extension, solver_types):
    source_file = output_dir + project.solver_name + s_ext + '.cpp'

    f = open(source_file, "r")
    lines = f.readlines()
    f.close()

    # Modify eigenvalues, fluxes, sources and initial conditions
    write_eigenvalues(lines, n_dust, c, s_type)
    write_flux(lines, n_dust, c)
    write_source(lines, n_dust, q, Stokes, weights, eta)
    write_initial(lines, n_dust, mu, Stokes, eta, s_type, Kx, Kz,
                  sigma=sigma)
    write_outflow_boundary(lines, project.n_vars, s_type)

    # Write to file
    f = open(source_file, "w")
    f.writelines(lines)
    f.close()

    # Write empty functions (no periodic boundaries by default)
    # Can only do this if periodic boundaries have been enabled
    if allow_periodic == True:
        write_periodic_dummies(output_dir, project.solver_name, s_ext)

# For a limiting scheme, need to specify whether DG solution is physical
if (project.solver_type == 'Limiting-ADER-DG'):
    write_physical(output_dir, project.solver_name, n_dust)


################################
# Implement periodic boundaries
################################
if args.periodic:
    print("Implementing periodic boundaries in solver...",
          project.solver_type)

    n_ghost = 1
    if project.solver_type == 'Limiting-ADER-DG':
        n_ghost = 2

    # Add Plot/Adjust periodic functions to solver class
    if (project.solver_type == 'ADER-DG'):
        write_periodic_functions(project.n_vars,
                                 project.order,
                                 [project.offset_x, project.offset_y],
                                 [project.size_x, project.size_y],
                                 output_dir,
                                 project.solver_name,
                                 n_ghost)
    elif (project.solver_type == 'Finite-Volumes'):
        write_periodic_functions(project.n_vars,
                                 project.patch_size,
                                 [project.offset_x, project.offset_y],
                                 [project.size_x, project.size_y],
                                 output_dir,
                                 project.solver_name,
                                 n_ghost)
    elif (project.solver_type == 'Limiting-ADER-DG'):
        # Modify the ADER-DG part, FV part not used
        write_periodic_functions(project.n_vars,
                                 project.order,
                                 [project.offset_x, project.offset_y],
                                 [project.size_x, project.size_y],
                                 output_dir,
                                 project.solver_name + '_ADERDG',
                                 n_ghost)
        write_periodic_functions(project.n_vars,
                                 2*project.order + 1,
                                 [project.offset_x, project.offset_y],
                                 [project.size_x, project.size_y],
                                 output_dir,
                                 project.solver_name + '_FV',
                                 n_ghost)
    else:
        print("Periodic boundary conditions for {} not implemented!".format(solver_type))
        exit(1)
