"""
Methods used across all or most workflows including:

    get_version: figure out the version of a command
    parse_stats: get the read and base counts from prinseq output
    apply_defaults: set dict (usually config) defaults recursively

"""

import os
import re
import logging
import subprocess
import tempfile
import yaml
import pandas
from snakemake.logging import logger

def get_version(command, version_flag='--version',
                cmd_prefix='',
                lines=None,
                regular_expression=None):
    """
    Gets the version string from a command

    cmd_prefix is useful if you need an interpreter (bash, python, etc)

    lines can be a line number (int), a slice object, or an iterable
    indicating which lines of the output to grab

    if regular_expression given, the first captured group is returned.
    """
    command = " ".join([cmd_prefix,
                        command,
                        version_flag,
                        "; exit 0"])
    logger.debug("Running command:\n{}".format(command))
    out = subprocess.check_output('bash -c "{}"'.format(command),
                                  stderr=subprocess.STDOUT,
                                  shell=True).decode()
    if re.search(r' not found', out):
        raise Exception("Command '{}' was not found:\n{}".format(command, out))


    # select specific lines
    if lines is not None:
        out_lines = out.split("\n")
        if isinstance(lines, slice):
            out = "\n".join(out_lines[lines])
        elif isinstance(lines, int):
            try:
                out = out_lines[lines]
            except IndexError:
                logger.warning(
                    "Line {} does not exist in output from '{}':\n{}"\
                                            .format(lines, command, out)
                )
        else:
            out = "\n".join(out_lines[i] for i in lines)

    # apply regular expression if given
    if regular_expression is None:
        return re.sub(r'[\n\r]',' ', out.strip())
    else:
        if isinstance(regular_expression, str):
            regular_expression = re.compile(regular_expression)
        match = regular_expression.search(out)
        if match is None:
            logger.warning(
                "Expression {} did not matach a group in output from `{}`: {}"\
                        .format(regular_expression.pattern, command, out)
            )
        else:
            try:
                return match.group(1)
            except IndexError:
                # expression mathed but no group specified, return entire match
                return match.group()

def parse_stats(stats_file):
    """
    pull out the read and base counts from a prinseq output file.
    Returns a two item dict with integer values and keys: 'reads', 'bases'
    """

    # if the file is empty, so was the fasta/fastq file
    if os.stat(stats_file).st_size == 0:
        return {'reads':0, 'bases':0}

    stats = pandas.read_table(stats_file,
                              names=('module', 'key', 'value'),
                              index_col=1)['value']
    return {k:int(stats[k]) for k in ['reads', 'bases']}

def add_stats_outputs(snakefile, config):
    """
    if config[discover_fastx_for_stats] is set to true,
    figure out what fasts[aq] files
    this workflow generates and add the T.stats and T.hist to the outputs for
    each of file T.

    We do this by calling the snakemake API to get the output summary.
    """
    if config.get('discover_fastx_for_stats', False) in [True, 'True']:
        config_copy = dict(config)
        config_copy['discover_fastx_for_stats'] = False
        config_file = tempfile.NamedTemporaryFile(mode='w')
        yaml.dump(config_copy, config_file)

        command = [
            'snakemake',
            '-s',
            snakefile,
            '--nolock',
            '--configfile',
            config_file.name,
            '--summary',
            '--rerun-incomplete',
            '-n',
        ]

        print("logger: {}\n {}".format(logger, dir(logger)))
        for item in dir(logger):
            print("logger.{}".format(item))
            try:
                print(eval("logger.{}()".format(item)))
            except:
                print(eval("logger.{}".format(item)))
        if logger.logger.getEffectiveLevel() >= logging.DEBUG:
            command += ['--verbose']

        logger.debug("Performing dry-run to get outputs")
        logger.debug(" ".join(command))
        try:
            complete = subprocess.run(command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        except:
            logger.warning("Cannot get fastx files, there is something wrong "
                            "with your workflow!")
            raise
            return

        if complete.returncode != 0:
            logger.warning("STDOUT:\n" + get_str(complete.stdout))
            logger.warning("STDERR:\n" + get_str(complete.stderr))
            raise Exception("Cannot get fastx files, there is something wrong "
                            "with your workflow!")

        logger.debug("Dry run complete")

        new_outputs = config.setdefault('outputs', set())
        output_count = 0
        for line in complete.stdout.decode().split('\n'):
            output = line.split('\t')[0]
            logger.debug(output)
            if re.search(r'f(aa|fn|na|a|asta|astq)(\.gz)?$', output):
                output_count += 1
                for extension in ['.hist', '.stats']:
                    new_outputs.add('stats/' + output + extension)

        logger.debug("Added stats and hist files for {} fasta files"\
                     .format(output_count))

def get_str(possibly_byte_array):
    if isinstance(possibly_byte_array, str):
        return possibly_byte_array
    return possibly_byte_array.decode()

def apply_defaults(config, defaults):
    """ recursively appy defaults to nested dicts """
    for param, pdefaults in defaults.items():
        if isinstance(pdefaults, dict):
            apply_defaults(config.setdefault(param, {}), pdefaults)
        else:
            config.setdefault(param, pdefaults)

def get_file_name(list_or_string):
    """
    Different versions of snakemake seem to have different approaches to this syntax:

    rule x:
        input:
            file1="some.file"
            file2="other.file"

    Sometimes input.file1 is a string and sometimes it is list containing a string.

    So I'm inserting this method into all the python code that uses named files.
    """
    if isinstance(list_or_string, str):
        return list_or_string
    return next(iter(list_or_string))
