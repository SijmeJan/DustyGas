from common import replace_with_indent

def write_eigenvalues(lines, n_dust, c):
    for i in range(0, len(lines)):
        # Gas normal velocity
        vg = 'Q[dIndex + 1]/Q[0]'

        # Gas eigenvalues
        lines[i] = replace_with_indent(lines[i],
                                      'lambda[0] = ',
                                      'lambda[0] = {} - {};\n'.format(vg, c))
        lines[i] = replace_with_indent(lines[i],
                                      'lambda[1] = ',
                                      'lambda[1] = {};\n'.format(vg))
        lines[i] = replace_with_indent(lines[i],
                                      'lambda[2] = ',
                                      'lambda[2] = {};\n'.format(vg))
        lines[i] = replace_with_indent(lines[i],
                                      'lambda[3] = ',
                                      'lambda[3] = {} + {};\n'.format(vg, c))

        # Dust eigenvalues
        for n in range(0, n_dust):
            for j in range(0, 4):
                l = 4*n + 4 + j

                # Dust normal velocity
                vd = 'Q[dIndex + {}]/Q[{}]'.format(5 + 4*n, 4 + 4*n)
                lines[i] =replace_with_indent(lines[i],
                                             'lambda[{}] = '.format(l),
                                             'lambda[{}] = {};\n'.format(l,vd))
