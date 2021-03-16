from common import replace_with_indent

def write_flux(lines, n_dust, c):
    for i in range(0, len(lines)):
        # First direction gas
        lines[i] = replace_with_indent(lines[i],
                                      'F[0][0] = ',
                                      'F[0][0] = Q[1];\n')
        f = 'F[0][1] = Q[1]*Q[1]/Q[0] + {}*Q[0];\n'.format(c*c)
        lines[i] = replace_with_indent(lines[i], 'F[0][1] = ', f)
        lines[i] = replace_with_indent(lines[i],
                                      'F[0][2] = ',
                                      'F[0][2] = Q[1]*Q[2]/Q[0];\n')
        lines[i] = replace_with_indent(lines[i],
                                      'F[0][3] = ',
                                      'F[0][3] = Q[1]*Q[3]/Q[0];\n')

        # Second direction gas
        lines[i] = replace_with_indent(lines[i],
                                      'F[1][0] = ',
                                      'F[1][0] = Q[2];\n')
        lines[i] = replace_with_indent(lines[i],
                                      'F[1][1] = ',
                                      'F[1][1] = Q[2]*Q[1]/Q[0];\n')
        f = 'F[1][2] = Q[2]*Q[2]/Q[0] + {}*Q[0];\n'.format(c*c)
        lines[i] = replace_with_indent(lines[i], 'F[1][2] = ', f)
        lines[i] = replace_with_indent(lines[i],
                                      'F[1][3] = ',
                                      'F[1][3] = Q[2]*Q[3]/Q[0];\n')

        # Dust fluxes
        for n in range(0, n_dust):
            # First direction
            lines[i] = replace_with_indent(lines[i],
                                           'F[0][{}] = '.format(4*n + 4),
                                           'F[0][{}] = Q[{}];\n'.format(4*n + 4, 5 + 4*n))
            lines[i] = replace_with_indent(lines[i],
                                           'F[0][{}] = '.format(4*n + 5),
                                           'F[0][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 5, 4*n + 5, 4*n + 5, 4*n + 4))
            lines[i] = replace_with_indent(lines[i],
                                           'F[0][{}] = '.format(4*n + 6),
                                           'F[0][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 6, 4*n + 5, 4*n + 6, 4*n + 4))
            lines[i] = replace_with_indent(lines[i],
                                           'F[0][{}] = '.format(4*n + 7),
                                           'F[0][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 7, 4*n + 5, 4*n + 7, 4*n + 4))

            # Second direction
            lines[i] = replace_with_indent(lines[i],
                                           'F[1][{}] = '.format(4*n + 4),
                                           'F[1][{}] = Q[{}];\n'.format(4*n + 4, 4*n + 6))
            lines[i] = replace_with_indent(lines[i],
                                           'F[1][{}] = '.format(4*n + 5),
                                           'F[1][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 5, 4*n + 6, 4*n + 5, 4*n + 4))
            lines[i] = replace_with_indent(lines[i],
                                           'F[1][{}] = '.format(4*n + 6),
                                           'F[1][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 6, 4*n + 6, 4*n + 6, 4*n + 4))
            lines[i] = replace_with_indent(lines[i],
                                           'F[1][{}] = '.format(4*n + 7),
                                           'F[1][{}] = Q[{}]*Q[{}]/Q[{}];\n'.format(4*n + 7, 4*n + 6, 4*n + 7, 4*n + 4))
