import subprocess
import glob
import re
import pandas


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


def parse_stats(stats_file):
    """
    pull out the read and base counts from a prinseq output file.
    Returns a two item dict with integer values and keys: 'reads', 'bases'
    """
    stats = pandas.read_table(stats_file,names=('module','key','value'),index_col=1)['value']
    return {k:int(stats[k]) for k in ['reads','bases']}


def collect_sample_reads(config):
    """
    If there is not already a "reads" dict in config, build it from the samples
    pattern
    """
    if 'reads' not in config:
        try:
            samples_pattern = config['samples_pattern']
        except:
            raise Exception("Please supply list of samples or patterns to find"
                            "them with") 
        reads = config.setdefault('reads',{})
        sample_RE = re.compile(samples_pattern.get('re',
                                                   r'/([^/]+)/[^/]+$'))
        #print(sample_RE.pattern)
        for read_file in glob.glob(samples_pattern.get('glob',
                                                       './*/reads.cleaned.fastq.gz')):
            # print(read_file)
            sample = sample_RE.search(read_file).group(1)
            # sanitize sample name
            sample = re.sub(r'[^A-Za-z0-9_]','_',sample)
            # print(sample)
            reads[sample]=read_file
#
