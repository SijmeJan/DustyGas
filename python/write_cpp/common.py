def replace_with_indent(line, search, replace):
    ret = line

    if (line.find(search) != -1):
        n_white  = len(line) - len(line.lstrip())
        ret = replace.rjust(len(replace) + n_white)

    return ret

def remove_function_body(lines, func_name):
    body = []
    for i in range(0, len(lines)):
        if (lines[i].find(func_name + '(') != -1):
            while (lines[i].find('{') == -1):
                i = i + 1

            # Deal with case of single line function definition: { return 0; }
            before_curly = lines[i].split("{",1)[0]
            after_curly = lines[i].split("{",1)[1]
            if len(after_curly) > 0:
                lines[i] = before_curly + '{\n'
                if after_curly.find("}") != -1:
                    after_curly = after_curly.split("}",1)[0] + '\n'
                    lines[i+1:i+1] = ['}\n']
                lines[i+1:i+1] = [after_curly]

            open_brackets = 1
            while (open_brackets > 0):
                if (lines[i+1].find('{') != -1):
                    open_brackets = open_brackets + 1
                if (lines[i+1].find('}') != -1):
                    open_brackets = open_brackets - 1
                if (open_brackets > 0):
                    body.append(lines[i+1])
                    lines.pop(i+1)
            break;

    return body

def add_function_body(lines, func_name, body):
    for i in range(0, len(lines)):
        if (lines[i].find(func_name + '(') != -1):
            while (lines[i].find('{') == -1):
                i = i + 1
            lines[i+1:i+1] = body
            break;

def add_class_member(lines, class_name, public, member):
    for i in range(0, len(lines)):
        if (lines[i].find(class_name) != -1 and
            lines[i].split()[0] == 'class'):
            while (lines[i].find(public) == -1):
                i = i + 1
            lines[i+1:i+1] = member
            break;
