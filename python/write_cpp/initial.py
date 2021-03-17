from common import replace_with_indent, remove_function_body

def write_initial(lines, n_dust, mu, Stokes):
    remove_function_body(lines, 'adjustSolution')

    tau = Stokes[0]
    denom = '((1 + {})*(1 + {}) + {}*{})'.format(mu, mu, tau, tau)

    initial = ['  if (t == 0.0) {\n',
               '    Q[0] = 1.0;\n',
               '    Q[1] = 2*{}*{}/{};\n'.format(mu, tau, denom),
               '    Q[2] = 0.0;\n',
               '    Q[3] = -(1 + {}*{}*{}/{})/(1 + {});\n'.format(mu, tau, tau, denom, mu)]
    for n in range(0, n_dust):
        initial.extend(['    Q[{}] = {};\n'.format(4*n + 4, mu),
                        '    Q[{}] = -2*{}*{}/{};\n'.format(4*n + 5, mu, tau, denom),
                        '    Q[{}] = 0.0;\n'.format(4*n + 6),
                        '    Q[{}] = -{}*(1-{}*{}/{})/(1+{});\n'.format(4*n + 7, mu, tau, tau, denom, mu)])
    initial.extend(['  }\n'])

    for i in range(0, len(lines)):
        if (lines[i].find('adjustSolution(') != -1):
            lines[i+1:i+1] = initial
            break;
