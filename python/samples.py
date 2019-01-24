"""
METHODS for setting up workflows with multiple samples

    collect_sample_reads: find read fastq files using glob/re
"""
import re
import os
import glob
import snakemake
from python.download import remote_wrapper

def process_sample_data(sample_data, config):
    """
    sample_data is a top level config map that has two types of entries

     - samples: keyed on sample name and containing paths and other data
     - 'reads_patterns': list of patterns to find samples. These are passed to 
            collect_sample_reads() below

    This method processes the read patterns to find samples

    AND

    returns a list of sample names
    """

    # First, process the patterns
    if 'reads_patterns' in sample_data:
        if isinstance(sample_data['reads_patterns'], dict):
            reads_patterns = [sample_data['reads_patterns']]
        else:
            reads_patterns = sample_data['reads_patterns']
        for pattern_data in reads_patterns:
            if pattern_data.get('cleaned', True) in [True, 'True']:
                read_key = 'clean'
            else:
                read_key = 'raw'
            for sample, reads in collect_sample_reads(pattern_data, config).items():
                sample_data.setdefault(sample, {})[read_key] = reads

        # Now get rid of any patterns from config
        del sample_data['reads_patterns']

    # return sample names
    return [s for s in sample_data if s != 'reads_patterns']


def collect_sample_reads(samples_pattern_data, config):
    """
    Use the samples_pattern entry in config to dynamically locate read files
    and group into samples

    There are two ways to specify a pattern. You can use snakemake style
    wildcard_glob:

    sample_data:
        samples_pattern:
            wildcard_glob: "../data/{sample}_{direction}.fastq"

    If the string begins with SFTP, the makefile will use the SFTP remote
    to get files from a remote system.

    Note: although snakemake generally allows wildcards to span directories,
    this implementation will ignore such wildcards.

    You can also use a bash style globto find files and regular expression to get sample names from them:

    sample_data:
        samples_pattern:
            glob: "../data/*.fastq"
            re: "/([^_]+)_[^/]+\\.fastq"

    The glob finds files using a
    filesystem wildcard and the re should identify (as the first matched group)
    the sample name in the found file names.

    Returns dict mapping from sample names to read files.

    In the above example, let the generated reads dict could look like:
    reads:
        sample-01: reads/sample-01/reads.corrected.bfc.fastq.gz
        sample-02: reads/sample-02/reads.corrected.bfc.fastq.gz

    or this (depending on filesystem contents):
    reads:
        sample-01:
            - reads/sample-01/reads.R1.fastq
            - reads/sample-01/reads.R2.fastq
        sample-02:
            - reads/sample-02/reads.R1.fastq
            - reads/sample-02/reads.R2.fastq

    You can specify multiple patterns and mix styles:

    sample_data:
        samples_pattern:
            - wildcard_glob: "../data/{sample}_{direction}.fastq"
            - glob: "../data/*.fastq"
              re: "/([^_]+)_[^/]+\\.fastq"

    """

    # setup
    # new way
    if 'wildcard_glob' in samples_pattern_data:
        return get_wc_glob_reads(samples_pattern_data['wildcard_glob'], config)

    # old way (deprecated)
    snakemake.logger.warning("getting reads from glob and re is deprecated, "
                   "use wildcard_glob instead!")
    sample_pattern = samples_pattern_data.get('re', r'/([^/]+)/[^/]+$')
    sample_RE = re.compile(sample_pattern)
    read_file_glob = samples_pattern_data.get('glob',
                                              './*/reads.cleaned.fastq.gz')
    return get_re_glob_reads(read_file_glob, sample_RE)

def get_wc_glob_reads(wildcard_glob_string, config):
    """
    Use wildcard glob string to build map from samples to read fastq files
    
    There must be a "sample" wildcard in the glob string.
    
    Other wildcards are ignored. Multiple fastq files per samples are assumed to be fwd/rev pair (in alphabetical order)
    """
    reads = {}

    wildcard_values = remote_wrapper(wildcard_glob_string, config, glob=True)
 
    # which wildcard in glob was "sample"
    wildcard_names = list(wildcard_values._fields)
    sample_index = wildcard_names.index('sample')

    # build reads dictionary
    for wildcard_tuple in zip(*list(wildcard_values)):
        for value in wildcard_tuple:
            # skip if wildcard spans folders
            if re.search(os.path.sep, value):
                break
        else:
            # add sample
            sample = wildcard_tuple[sample_index]
            reads.setdefault(sample, []).append(
                remote_wrapper(wildcard_glob_string \
                               .format(**dict(zip(wildcard_names,
                                                  wildcard_tuple))),
                               config))
    return reads

def get_re_glob_reads(read_file_glob, sample_RE):
    """
    Given a glob string and regular expression,
    find all the matching files and use rexp
    to get a sample name for each file.

    return a dictionary from samples to lists of files
    """
    reads = {}
    read_files = glob.glob(read_file_glob)
    if len(read_files) == 0:
        raise Exception(
            "The sample reads wildcard '{}' did not match any files!"\
                            .format(read_file_glob)
        )

    # collect files into lists by sample
    for read_file in read_files:
        match = sample_RE.search(read_file)
        if match is None:
            raise Exception(
                ("The sample matching expression ({}) failed to find a sample "
                 "name in the path: {}").format(sample_pattern, read_file)
            )
        sample = match.group(1)
        # sanitize sample name
        sample = re.sub(r'[^A-Za-z0-9_]', '_', sample)
        reads.setdefault(sample, []).append(read_file)

    return reads

def get_sample_reads_for_mapping(wildcards, config):
    """
    reads for mapping can come from different places
    
    if a for_mapping file is specified for this sample, use it
    if map_clean_reads is set to True, use the clean reads

    if nothing explicit is specified:
    use raw reads preferentially
     (using the renamed.interleaved.fastq version if we can)
    fall back to clean reads
    """
    if isinstance(wildcards, str):
        # allow for a diret call with the sample name
        sample = wildcards
    else:
        # normally just get sample name from wildcards object
        sample = wildcards.sample

    # look in configuration first
    if sample in config['sample_data'] and (
        'raw' in config['sample_data'][sample] or \
        'for_mapping' in config['sample_data'][sample] or \
        'clean' in config['sample_data'][sample]
    ):

        if 'for_mapping' in config['sample_data'][sample]:
            # use explicit decalration first
            return config['sample_data'][sample]['for_mapping']

        if config.get('map_clean_reads', False):
            # return clean reads if explicitly requested
            return config['sample_data'][sample]['clean']

        try:
            # can we find the renamed and interleaved reads?
            renamed_root = re.search(r'^.+\.renamed(?:\.interleaved)?',
                                     config['sample_data']\
                                                    [sample]\
                                                    ['clean']
                                    ).group()
            return renamed_root + ".fastq"
        except:
            # OK, that didn't work
            pass
        # let's try simply using raw or clean reads
        if 'raw' in config['sample_data'][sample]:
            # just use the raw reads
            return config['sample_data'][sample]['raw']
        else:
            # return clean reads if that's all we have
            return config['sample_data'][sample]['clean']
    else:
        # otherwise, look for fastq file with sample name
        return "{}.fastq".format(sample)


