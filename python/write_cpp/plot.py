from common import replace_with_indent, remove_function_body, add_function_body

def write_plotter_gas_velocity(lines):
    remove_function_body(lines, 'mapQuantities')

    body = ['  const int writtenUnknowns = 3;\n',
            '  for (int i=0; i<writtenUnknowns; i++){\n',
            '    outputQuantities[i] = Q[i+1];\n',
            '  }\n']

    add_function_body(lines, 'mapQuantities', body)
