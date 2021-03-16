from common import replace_with_indent, remove_function_body

def write_initial(lines, n_dust, q, Stokes):
    remove_function_body(lines, 'adjustSolution')

    initial = ['  if (t == 0.0) {\n',
               '    Q[0] = 1.0;\n',
               '    Q[1] = 0.0;\n',
               '    Q[2] = 0.0;\n',
               '    Q[3] = 0.0;\n',
               '    if (sqrt((x[0] - 0.7)*(x[0] - 0.7) + (x[1] - 0.7)*(x[1] - 0.7)) < 0.2) Q[0] = 2.0;\n',
               '  }\n']

    for i in range(0, len(lines)):
        if (lines[i].find('adjustSolution(') != -1):
            lines[i+1:i+1] = initial
            break;
