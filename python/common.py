import subprocess


def get_version(command, version_flag='--version', 
                cmd_prefix='',
                lines=None,
                regular_expression=None):
    """
    Gets the version string from a command

    cmd_prefix is useful if you need an interpreter (bash, python, etc)
    
    lines can be a line number (int), a slice object, or an iterable indicating which lines of the output to grab

    if regular_expression given, the first captured group is returned.
    """
    out = subprocess.check_output(" ".join([cmd_prefix,
                                            command,
                                            version_flag,
                                            "; exit 0"]),
                                  stderr=subprocess.STDOUT,
                                  shell=True).decode()

    # select specific lines
    if lines is not None:
        out_lines = out.split("\n")
        if isinstance(lines,slice):
            out = "\n".join(out_lines[lines])
        elif isinstance(lines, int):
            out = out_lines[lines]
        else:
            out = "\n".join(out_lines[i] for i in lines)

    # apply regular expression if given
    if regular_expression is None:
        return out.strip()
    else:
        return regular_expression.search(out).group(1)
