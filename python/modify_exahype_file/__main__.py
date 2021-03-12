import argparse

class LoadFromFile (argparse.Action):
    def __call__ (self, parser, namespace, values, option_string = None):
        with values as f:
            # parse arguments in the file and store them in the target namespace
            parser.parse_args(f.read().split(), namespace)

def replace_with_indent(line, search, replace):
    ret = line

    if (line.find(search) != -1):
        n_white  = len(line) - len(line.lstrip())
        ret = replace.rjust(len(replace) + n_white)

    return ret

parser = argparse.ArgumentParser(description='Create a DustyGas simulation')

parser.add_argument('infile',
                    help='ExaHyPE file to modify')
parser.add_argument('--project-name', dest='project_name', type=str,
                    default='MyProject',
                    help='Name of the ExaHyPE project')
parser.add_argument('--n_dust', dest='n_dust', type=int,
                    default=0,
                    help='Number of dust fluids')
parser.add_argument('--dims', dest='dims', type=int,
                    default=2,
                    help='Number of space dimensions')
parser.add_argument('--box-size', dest='bsize', type=float, nargs='+',
                    default=[1.0, 1.0],
                    help='Box size x and y')
parser.add_argument('--outfile', dest='outfile', type=str,
                    default='dustygas.exahype',
                    help='Output exahype file')
parser.add_argument('--parfile', type=open, action=LoadFromFile,
                    help='Parameter file')
args = parser.parse_args()

if (len(args.bsize) != args.dims):
    print('Error: Number of elements in domain size must equal number of spatial dimensions ({}), but got {}'.format(args.dims, args.bsize))
    exit(1)

# Total number of unknowns; gas + dust
n_unknowns = 4 + 4*args.n_dust

# Bottom left corner of domain
xL = -0.5*args.bsize[0]
yL = -0.5*args.bsize[1]

f = open(args.infile, "r")
lines = f.readlines()
f.close()

for i in range(0, len(lines)):
    lines[i] = replace_with_indent(lines[i],
                                   'exahype-project ',
                                   'exahype-project ' +
                                       args.project_name + '\n')
    lines[i] = replace_with_indent(lines[i],
                                   'variables const',
                                   'variables const = {}\n'.format(n_unknowns))
    lines[i] = replace_with_indent(lines[i],
                                   'dimension const',
                                   'dimension const = {}\n'.format(args.dims))
    lines[i] = replace_with_indent(lines[i],
                                   'width',
                                   'width = {}, {}\n'.format(args.bsize[0],
                                                             args.bsize[1]))
    lines[i] = replace_with_indent(lines[i],
                                   'offset',
                                   'offset = {}, {}\n'.format(xL, yL))

f = open(args.outfile, "w")
f.writelines(lines)
f.close()
