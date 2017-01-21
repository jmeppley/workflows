import subprocess
import glob
import re
import pandas
import os


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
    command = " ".join([cmd_prefix,
                        command,
                        version_flag,
                        "; exit 0"])
    out = subprocess.check_output(command,
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
        if isinstance(regular_expression, str):
            regular_expression = re.compile(regular_expression)
        try:
            return regular_expression.search(out).group(1)
        except AttributeError:
            print("WARNING: Expression {} did not matach a group in output from `{}`: {}".format(regular_expression.pattern, command, out))


def parse_stats(stats_file):
    """
    pull out the read and base counts from a prinseq output file.
    Returns a two item dict with integer values and keys: 'reads', 'bases'
    """

    # if the file is empty, so was the fasta/fastq file
    if os.stat(stats_file).st_size==0:
        return {'reads':0,'bases':0}

    stats = pandas.read_table(stats_file,names=('module','key','value'),index_col=1)['value']
    return {k:int(stats[k]) for k in ['reads','bases']}

def apply_defaults(config, defaults):
    """ recursively appy defaults to nested dicts """
    for param, pdefaults in defaults.items():
        if isinstance(pdefaults, dict):
            apply_defaults(config.setdefault(param,{}), pdefaults)
        else:
            config.setdefault(param, pdefaults)

def collect_sample_reads(config, get_stats=True):
    """
    If there is not already a "reads" dict in config, build it from the samples
    pattern

    returns True if it sets up filenames that need QC rules.

    the resulting reads dict should look like this:
    config['reads']={
        'sample_1': '/path/to/sample_1.fastq',
        'sample_2': '/path/to/sample_2.fastq'
    }
    
    The indicated samples files may exist somewhere else or may be targets that 
    need to be built. In the latter case, there should be transitions set to
    define how the base files is to be generated. For example, we are starting 
    with raw reads, the reads above will be the cleaned versions and snakemake
    will have to figure out how to generate them.

    The samples_pattern dict should look like this:
        samples_pattern:
            glob: "../data/*.fastq"
            re: "/([^_]+)_[^/]+\.fastq"
            cleaned: False
            name: all
            metadata: samples.all.tsv

    In the above example, let the generated reads dict could look like:
    reads:
        sample-01: reads/sample-01/reads.corrected.bfc.fastq.gz
        sample-02: reads/sample-02/reads.corrected.bfc.fastq.gz

    Snakemake understands that reads/{sample}/reads.corrected.bfc.fastq.gz
    is generated from a series of steps from reads/{sample}/reads.R1.fastq and
    reads/{sample}/reads.R2.fastq. Those files don't exist, so 4 transitions
    will also be defined:
    transitions:
        reads/sample-01/reads.R1.fastq: ../data/sample-01_R1.fastq
        reads/sample-01/reads.R2.fastq: ../data/sample-01_R2.fastq
        reads/sample-02/reads.R1.fastq: ../data/sample-02_R1.fastq
        reads/sample-02/reads.R2.fastq: ../data/sample-02_R2.fastq
    """

    # don't do anything if reads are already set
    if 'reads' in config:
        return False

    # initialize and check some basics
    try:
        samples_pattern_data = config['samples_pattern']
    except:
        raise Exception("Please supply list of samples or patterns to find"
                        "them with") 
    reads = config.setdefault('reads',{})
    sample_pattern = samples_pattern_data.get('re', r'/([^/]+)/[^/]+$')
    sample_RE = re.compile(sample_pattern)
    read_file_glob = samples_pattern_data.get('glob',
                                              './*/reads.cleaned.fastq.gz')

    # find files
    read_files = glob.glob(read_file_glob)
    if len(read_files)==0:
        raise Exception("The sample reads wildcard '{}' did not match any files!"\
                            .format(read_file_glob))

    # collect files into lists by sample
    for read_file in read_files:
        m = sample_RE.search(read_file)
        if m is None:
            raise Exception("The sample matching expression ({}) failed to find a sample name in the path: {}".format(sample_pattern, read_file))
        sample = m.group(1)
        # sanitize sample name
        sample = re.sub(r'[^A-Za-z0-9_]','_',sample)
        reads.setdefault(sample,[]).append(read_file)

    # loop back over samples and set up cleaning or interleaving if needed
    needs_qc = False
    transitions = config.setdefault('transitions',{})
    already_cleaned = samples_pattern_data.get('cleaned',True) \
                        in [True, 1, "True", "T", "true", "t"]
    fasta_files = []
    for sample in list(reads.keys()):
        files = sorted(reads[sample])

        # Bail out if we have too many files per sample
        if len(files)>2:
            raise Exception("I don't know how to deal with more than two"
                            " files per sample!\nSample={}\nFiles:\n{}"\
                                .format(sample,
                                        "\n".join(files)))

        # Do we need to clean?
        if not already_cleaned:
            needs_qc = True
            # set target that will invoke QC
            reads[sample] = \
                    'reads/{sample}/reads.renamed.R12.cleaned.corrected.fastq.gz'\
                            .format(sample=sample)

            # keep track of fasta files that will be generated
            fasta_files.extend([
               'reads/{sample}/reads.renamed.R12.cleaned.corrected.fastq.gz'\
                                                            .format(sample=sample),
               'reads/{sample}/reads.renamed.R12.cleaned.fastq.gz'\
                                                            .format(sample=sample)])
                
        else:
            # we just need to get interleaved reads
            reads[sample] = 'reads/{sample}/reads.R12.fastq'\
                            .format(sample=sample)

        # both tracks go through this file
        fasta_files.append('reads/{sample}/reads.renamed.R12.fastq'\
                                                            .format(sample=sample))
            
        # Do we have a pair of files or single file?
        if len(files)==1:
            # we are starting from interleaved, just link it in
            transitions['reads/{sample}/reads.renamed.R12.fastq'\
                                .format(sample=sample)] = files[0]
        else:
            needs_qc = True
            # we are starting from paired files, link those
            transitions['reads/{sample}/reads.R1.fastq'\
                            .format(sample=sample)] = files[0]
            transitions['reads/{sample}/reads.R2.fastq'\
                            .format(sample=sample)] = files[1]
            
    if get_stats:
        # add stats file to final outputs for each fasta file generated
        outputs = config.setdefault('outputs', set())
        for ext in ['stats', 'hist']:
            for file_name in fasta_files:
                outputs.add(".".join((file_name, ext)))

    return needs_qc
