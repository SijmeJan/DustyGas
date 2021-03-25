from common import replace_with_indent

def write_source(lines, n_dust, q, Stokes, eta):

    for i in range(0, len(lines)):
        # Source x: eta + Coriolis
        lines[i] = replace_with_indent(lines[i],
                                      'S[1] = ',
                                      'S[1] = 2*Q[0]*{} + 0*Q[3];\n'.format(eta))
        # Source y: Coriolis
        lines[i] = replace_with_indent(lines[i],
                                      'S[3] = ',
                                      'S[3] = 0*{}*Q[1];\n'.format(q - 2))

        for n in range(0, n_dust):
            # Source x: Coriolis + drag
            lines[i] = replace_with_indent(lines[i],
                                           'S[{}] = '.format(4*n + 5),
                                           'S[{}] = 2*Q[{}] - (Q[{}] - Q[{}]*Q[1]/Q[0])/{};\n'.format(4*n + 5, 4*n + 7, 4*n + 5, 4*n + 4, Stokes[n]))

            # Source z: drag
            lines[i] = replace_with_indent(lines[i],
                                           'S[{}] = '.format(4*n + 6),
                                           'S[{}] = - (Q[{}] - Q[{}]*Q[2]/Q[0])/{};\n'.format(4*n + 6, 4*n + 6, 4*n + 4, Stokes[n]))

            # Source y: Coriolis + drag
            lines[i] = replace_with_indent(lines[i],
                                          'S[{}] = '.format(4*n + 7),
                                          'S[{}] = {}*Q[{}] - (Q[{}] - Q[{}]*Q[3]/Q[0])/{};\n'.format(4*n + 7, q - 2, 4*n + 5, 4*n + 7, 4*n + 4, Stokes[n]))

    gas_x_drag = []
    gas_y_drag = []
    gas_z_drag = []
    for n in range(0, n_dust):
        gas_x_drag.append('  S[1] += (Q[{}] -Q[{}]*Q[1]/Q[0])/{};\n'.format(4*n + 5, 4*n + 4, Stokes[n]))
        gas_z_drag.append('  S[2] += (Q[{}] -Q[{}]*Q[2]/Q[0])/{};\n'.format(4*n + 6, 4*n + 4, Stokes[n]))
        gas_y_drag.append('  S[3] += (Q[{}] -Q[{}]*Q[3]/Q[0])/{};\n'.format(4*n + 7, 4*n + 4, Stokes[n]))

    for i in range(0, len(lines)):
        if (lines[i].find('S[1] = ') != -1):
            lines[i+1:i+1] = gas_x_drag
    for i in range(0, len(lines)):
        if (lines[i].find('S[2] = ') != -1):
            lines[i+1:i+1] = gas_z_drag
    for i in range(0, len(lines)):
        if (lines[i].find('S[3] = ') != -1):
            lines[i+1:i+1] = gas_y_drag
