import numpy as np

from common import replace_with_indent, remove_function_body

class InitialConditions():
    def __init__(self):
        self.initial = ['  // To be implemented\n']

    def write(self, lines, solver_type):
        function_name = 'adjustSolution'
        if (solver_type == 'ADER-DG'):
            function_name = 'adjustPointSolution'

        remove_function_body(lines, function_name)

        for i in range(0, len(lines)):
            if (lines[i].find(function_name + '(') != -1):
                lines[i+1:i+1] = self.initial
                break;

class GasExplosionIC(InitialConditions):
    def __init__(self, centre, size, amp):
        self.initial = \
          ['  if (t == 0.0) {\n',
           '    Q[0] = 1.0;\n',
           '    Q[1] = 0.0;\n',
           '    Q[2] = 0.0;\n',
           '    Q[3] = 0.0;\n',
           '    if (sqrt((x[0] - {})**2 + (x[1] - {})**2) < {})\n'.format(centre[0], centre[1], size),
           '      Q[0] = {};\n'.format(1 + amp),
           '  }\n']

class GasDensityWaveIC(InitialConditions):
    def __init__(self, Kx, amp):
        omega = np.sqrt(Kx*Kx + 1.0)

        self.initial = \
          ['  if (t == 0.0) {\n',
           '    double a = {};\n'.format(amp),
           '    double c = cos(x[0]*{});\n'.format(Kx),
           '    double s = sin(x[0]*{});\n'.format(Kx),
           '\n',
           '    Q[0] = 1.0 + {}*a*c/{};\n'.format(Kx, omega),
           '    Q[1] = a*c;\n',
           '    Q[2] = 0.0;\n',
           '    Q[3] = 0.5*a*s/{};\n'.format(omega),
           '  }\n']


class MonoDustyGasIC(InitialConditions):
    def __init__(self, Kx, Kz, amp, mu, Stokes, n_dust, eta):
        if n_dust < 1:
            raise RuntimeError('Cannot initiate Dusty Gas with less than one dust fluid')
        tau = Stokes[0]
        denom = '((1 + {})*(1 + {}) + {}*{})'.format(mu, mu, tau, tau)

        initial = \
          ['  if (t == 0.0) {\n',
           '    double a = {};\n'.format(amp),
           '    double c = cos(x[0]*{} + x[1]*{});\n'.format(Kx, Kz),
           '    double s = sin(x[0]*{} + x[1]*{});\n'.format(Kx, Kz),
           '\n',
           '    Q[0] = 1.0;\n',
           '    Q[1] = 2*{}*{}/{};\n'.format(mu, tau, denom),
           '    Q[2] = 0.0;\n',
           '    Q[3] = -(1 + {}*{}*{}/{})/(1 + {});\n'.format(mu, tau, tau, denom, mu)]

        # Dust density and *velocities*
        for n in range(0, n_dust):
            initial.extend(\
                ['    Q[{}] = {};\n'.format(4*n + 4, mu),
                 '    Q[{}] = -2*{}/{};\n'.format(4*n + 5, tau, denom),
                 '    Q[{}] = 0.0;\n'.format(4*n + 6),
                 '    Q[{}] = -(1-{}*{}/{})/(1+{});\n'.format(4*n + 7, tau, tau, denom, mu)])

        # Add eigenvector (linearA)
        eigen = [+0.0000074637 + 0.0000070677*1j,
                 -0.0563787907 + 0.0120535455*1j,
                 +0.0563784989 - 0.0120536242*1j,
                 +0.0445570113 + 0.0197224299*1j,
                 1.0,
                 -0.0466198076 + 0.0124333223*1j,
                 +0.0546507401 - 0.0077776652*1j,
                 +0.0435211557 + 0.0213517453*1j]

        for n in range(0, 4 + 4*n_dust):
            initial.extend(['    Q[{}] += a*({}*c + {}*s);\n'.format(n, np.real(eigen[n]), -np.imag(eigen[n]))])

        # Change into momenta
        initial.extend(['    Q[1] *= Q[0];\n',
                        '    Q[2] *= Q[0];\n',
                        '    Q[3] *= Q[0];\n'])
        for n in range(0, n_dust):
            initial.extend(['    Q[{}] *= Q[{}];\n'.format(4*n + 5, 4*n + 4),
                            '    Q[{}] *= Q[{}];\n'.format(4*n + 6, 4*n + 4),
                            '    Q[{}] *= Q[{}];\n'.format(4*n + 7, 4*n + 4)])

        if (eta != 1.0):
            # Adjust velocities in terms of eta
            initial.extend(['    Q[1] *= {};\n'.format(eta),
                            '    Q[2] *= {};\n'.format(eta),
                            '    Q[3] *= {};\n'.format(eta)])
        for n in range(0, n_dust):
            initial.extend(['    Q[{}] *= {};\n'.format(4*n + 5, eta),
                            '    Q[{}] *= {};\n'.format(4*n + 6, eta),
                            '    Q[{}] *= {};\n'.format(4*n + 7, eta)])

        initial.extend(['  }\n'])

        self.initial = initial

class PolyDustyGasIC(InitialConditions):
    def __init__(self, Kx, Kz, amp, mu, Stokes, n_dust, eta, sigma):
        if n_dust < 2:
            raise RuntimeError('Cannot initiate polydispers Dusty Gas with less than two dust fluids')

        J0 = sigma.Jint(0, mu)
        J1 = sigma.Jint(1, mu)

        denom = (1 + J0)*(1 + J0) + J1*J1

        initial = \
          ['  if (t == 0.0) {\n',
           '    double a = {};\n'.format(amp),
           '    double c = cos(x[0]*{} + x[1]*{});\n'.format(Kx, Kz),
           '    double s = sin(x[0]*{} + x[1]*{});\n'.format(Kx, Kz),
           '\n',
           '    // Gas\n',
           '    Q[0] = 1.0;\n',
           '    Q[1] = {};\n'.format(2*J1/denom),
           '    Q[2] = 0.0;\n',
           '    Q[3] = {};\n'.format(-(1 + J0)/denom)]

        # Dust density and *velocities*
        for n in range(0, n_dust):
            initial.extend(\
                ['    // Dust, Stokes = {}\n'.format(Stokes[n]),
                 '    Q[{}] = {};\n'.format(4*n + 4, mu*sigma.sigma(Stokes[n])),
                 '    Q[{}] = {};\n'.format(4*n + 5, 2*(J1 - Stokes[n]*(1 + J0))/denom/(1 + Stokes[n]*Stokes[n])),
                 '    Q[{}] = 0.0;\n'.format(4*n + 6),
                 '    Q[{}] = {};\n'.format(4*n + 7, -(1 + J0 + Stokes[n]*J1)/denom/(1 + Stokes[n]*Stokes[n]))])

        # Add perturbations (random)
        initial.extend('    // Add perturbations\n')
        rng = np.random.default_rng(12345)
        for n in range(0, 4 + 4*n_dust):
            initial.extend(['    Q[{}] += a*({}*c + {}*s);\n'.format(n, rng.random() - 0.5, rng.random() - 0.5)])

        # Change into momenta
        initial.extend(['    // Change to momenta\n',
                        '    Q[1] *= Q[0];\n',
                        '    Q[2] *= Q[0];\n',
                        '    Q[3] *= Q[0];\n'])
        for n in range(0, n_dust):
            initial.extend(['    Q[{}] *= Q[{}];\n'.format(4*n + 5, 4*n + 4),
                            '    Q[{}] *= Q[{}];\n'.format(4*n + 6, 4*n + 4),
                            '    Q[{}] *= Q[{}];\n'.format(4*n + 7, 4*n + 4)])

        if (eta != 1.0):
            # Adjust velocities in terms of eta
            initial.extend(['    // Velocities in units of eta\n',
                            '    Q[1] *= {};\n'.format(eta),
                            '    Q[2] *= {};\n'.format(eta),
                            '    Q[3] *= {};\n'.format(eta)])
            for n in range(0, n_dust):
                initial.extend(['    Q[{}] *= {};\n'.format(4*n + 5, eta),
                                '    Q[{}] *= {};\n'.format(4*n + 6, eta),
                                '    Q[{}] *= {};\n'.format(4*n + 7, eta)])

        initial.extend(['  }\n'])

        self.initial = initial

def write_initial(lines, n_dust, mu, Stokes, eta, solver_type,
                  Kx, Kz, sigma=None):
    # Epicyclic oscillation
    #ic = GasDensityWaveIC(Kx=0.0, amp=0.1)

    # 1D gas density wave
    #ic = GasDensityWaveIC(Kx=Kx, amp=0.001)

    ic = GasExplosionIC([np.pi/Kx, np.pi/Kz], 0.25*2*np.pi/Kx, 1.0)

    # LinearA test
    #ic = MonoDustyGasIC(Kx=Kx,
    #                    Kz=Kz,
    #                    amp=0.001,
    #                    mu=mu,
    #                    Stokes=Stokes,
    #                    n_dust=n_dust,
    #                    eta=eta)

    # Polydisperse
    #ic = PolyDustyGasIC(Kx=Kx,
    #                    Kz=Kz,
    #                    amp=0.001,
    #                    mu=mu,
    #                    Stokes=Stokes,
    #                    n_dust=n_dust,
    #                    eta=eta,
    #                    sigma=sigma)

    ic.write(lines, solver_type)
