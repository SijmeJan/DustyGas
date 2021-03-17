from common import replace_with_indent, remove_function_body, add_function_body

def write_plotter_gas_velocity(lines):
    remove_function_body(lines, 'mapQuantities')

    body = ['  const int writtenUnknowns = 3;\n',
            '  for (int i=1; i<writtenUnknowns+1; i++){\n',
            '    outputQuantities[i] = Q[i];\n',
            '  }\n']

    add_function_body(lines, 'mapQuantities', body)
