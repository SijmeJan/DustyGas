from common import remove_function_body, add_function_body

def write_physical(output_dir, solver_name, n_dust):
    source_file = output_dir + 'Abstract' + solver_name + '_ADERDG.h'

    f = open(source_file, "r")
    lines = f.readlines()
    f.close()

    remove_function_body(lines, 'isPhysicallyAdmissible')

    body = ['  // Check if all quantities are finite\n']

    for i in range(0, 4 + 4*n_dust):
        body.append('  if (!std::isfinite(localDMPObservablesMin[{}])) return false;\n'.format(i))
        body.append('  if (!std::isfinite(localDMPObservablesMax[{}])) return false;\n'.format(i))

    body.append('\n  // Check for positive densities\n')
    for i in range(0, n_dust + 1):
        body.append('  if (localDMPObservablesMin[{}] <= 0.0) return false;\n'.format(4*i))

    body.append('  return true;\n')

    add_function_body(lines, 'isPhysicallyAdmissible', body)

    # Write to file
    f = open(source_file, "w")
    f.writelines(lines)
    f.close()
