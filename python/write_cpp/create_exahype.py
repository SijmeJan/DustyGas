import numpy as np
import argparse
import os
import configparser

from scipy.optimize import fsolve

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


# Needs single argument; INI file
parser = argparse.ArgumentParser(description='Create ExaHyPE file from DustyGas INI')
parser.add_argument('--periodic', action="store_true", help="use periodic domain")
parser.add_argument('infile', help='INI file')
args = parser.parse_args()

# Output file: same name, different extension
output_file = args.infile.split('.')[0] + '.exahype'

# Parse INI file
config = configparser.ConfigParser()
config.read(args.infile)

Kx = float(config['DOMAIN']['Kx'])
Kz = float(config['DOMAIN']['Kz'])

mesh_size = float(config['DOMAIN']['mesh_size'])
end_time = float(config['DOMAIN']['end_time'])

direc = config['OUTPUT']['direc']
save_dt = float(config['OUTPUT']['save_dt'])

solver_type = config['SOLVER']['type']
order = int(config['SOLVER']['order'])
mpi = config['SOLVER']['mpi']
if solver_type == 'Limiting-ADER-DG':
    limiter = config['SOLVER']['limiter']

n_dust = int(config['DUST']['n_dust'])

n_vars = 4 + 4*n_dust
width_x = 2*np.pi/Kx
width_y = 2*np.pi/Kz

# Create periodic domain by addding ghost cells
if args.periodic:
    # Number of ghost cells required
    n_ghost = 1
    if solver_type == 'Limiting-ADER-DG':
        n_ghost = 2

    domain_size = \
      periodic_domain_size([width_x, width_y], n_ghost, mesh_size)
    width_x = domain_size[0]
    width_y = domain_size[1]

exahype_body = [\
  'exahype-project DustyGas\n',
  '  peano-kernel-path const = ./Peano\n',
  '  exahype-path const = ./ExaHyPE\n',
  '  output-directory const = {}\n'.format(direc),
  '\n',
  '  computational-domain\n',
  '    dimension const = 2\n',
  '    width = {}, {}\n'.format(width_x, width_y),
  '    offset = 0.0, 0.0\n',
  '    end-time = {}\n'.format(end_time),
  '  end computational-domain\n']

if mpi == 'yes':
  exahype_body[-1:-1] = [\
    '  \n',
    '  distributed-memory\n',
    '    identifier               = static_load_balancing\n',
    '    configure                = {greedy-regular,FCFS,ranks-per-node:1}\n',
    '    buffer-size              = 64\n',
    '    timeout                  = 120\n',
    '  end distributed-memory\n']

exahype_body[-1:-1] = [\
  '\n',
  '  solver {} DustyGasSolver\n'.format(solver_type),
  '    variables const = {}\n'.format(n_vars),
  '    order const = {}\n'.format(order),
  '    maximum-mesh-size = {}\n'.format(mesh_size),
  '    time-stepping = global\n',
  '    type const = nonlinear\n',
  '    terms const = flux, source\n',
  '    optimisation const = generic\n',
  '    language const = C\n']

if solver_type == 'Limiting-ADER-DG':
  exahype_body[-1:-1] = [\
    '\n',
    '    limiter-type const = {}\n'.format(limiter),
    '    limiter-optimisation const = generic\n',
    '    limiter-language const = C\n',
    '    dmp-observables = {}\n'.format(n_vars),
    '    dmp-relaxation-parameter = 1e14\n',
    '    dmp-difference-scaling = 1e13\n']

exahype_body[-1:-1] = [\
  '\n',
  '    plot vtk::Cartesian::cells::ascii State\n',
  '      variables const = {}\n'.format(n_vars),
  '      time = 0.0\n',
  '      repeat = {}\n'.format(save_dt),
  '      output = ./state\n',
  '    end plot\n',
  '\n',
  '  end solver\n',
  'end exahype-project']

# Write to file
f = open(output_file, "w")
f.writelines(exahype_body)
f.close()
