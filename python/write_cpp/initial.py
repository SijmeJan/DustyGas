import numpy as np
from common import replace_with_indent, remove_function_body

def write_initial(lines, n_dust, mu, Stokes):
    remove_function_body(lines, 'adjustSolution')

    Kx = 30.0
    Kz = 30.0

    tau = Stokes[0]
    denom = '((1 + {})*(1 + {}) + {}*{})'.format(mu, mu, tau, tau)

    initial = ['  if (t == 0.0) {\n',
               '    double a = 0.001;\n',
               '    double c = cos(x[0]*{} + x[1]*{});\n'.format(Kx, Kz),
               '    double s = sin(x[0]*{} + x[1]*{});\n'.format(Kx, Kz),
               '\n',
               '    Q[0] = 1.0;\n',
               '    Q[1] = 2*{}*{}/{};\n'.format(mu, tau, denom),
               '    Q[2] = 0.0;\n',
               '    Q[3] = -(1 + {}*{}*{}/{})/(1 + {});\n'.format(mu, tau, tau, denom, mu)]

    # Dust density and *velocities*
    for n in range(0, n_dust):
        initial.extend(['    Q[{}] = {};\n'.format(4*n + 4, mu),
                        '    Q[{}] = -2*{}/{};\n'.format(4*n + 5, tau, denom),
                        '    Q[{}] = 0.0;\n'.format(4*n + 6),
                        '    Q[{}] = -(1-{}*{}/{})/(1+{});\n'.format(4*n + 7, tau, tau, denom, mu)])

    # Add eigenvector
    eigen = [0.0000224 + 0.0000212*1j,
             -0.1691398 + 0.0361553*1j,
             0.1691389 - 0.0361555*1j,
             0.1336704 + 0.0591695*1j,
             1.0,
             -0.1398623 + 0.0372951*1j,
             0.1639549 - 0.0233277*1j,
             0.1305628 + 0.0640574*1j]

    for n in range(0, 4 + 4*n_dust):
        initial.extend(['    Q[{}] += a*({}*c + {}*s);\n'.format(n, np.real(eigen[n]), -np.imag(eigen[n]))])

    # Change into momenta
    initial.extend = ['    Q[1] *= Q[0];\n',
                      '    Q[2] *= Q[0];\n',
                      '    Q[3] *= Q[0];\n']
    for n in range(0, n_dust):
        initial.extend(['    Q[{}] *= Q[{}];\n'.format(4*n + 5, 4*n + 4),
                        '    Q[{}] *= Q[{}];\n'.format(4*n + 6, 4*n + 4),
                        '    Q[{}] *= Q[{}];\n'.format(4*n + 7, 4*n + 4)])

    initial.extend(['  }\n'])


    for i in range(0, len(lines)):
        if (lines[i].find('adjustSolution(') != -1):
            lines[i+1:i+1] = initial
            break;
