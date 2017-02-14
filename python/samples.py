"""
METHODS for setting up workflows with multiple samples

    collect_sample_reads: find read fastq files using glob/re
"""
import re
import glob

def collect_sample_reads(samples_pattern_data):
    """
    Use the samples_pattern entry in config to locate read files and group into
    samples

    The samples_pattern dict should look like this:
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

    """

    # setup
    reads = {}
    sample_pattern = samples_pattern_data.get('re', r'/([^/]+)/[^/]+$')
    sample_RE = re.compile(sample_pattern)
    read_file_glob = samples_pattern_data.get('glob',
                                              './*/reads.cleaned.fastq.gz')

    # find files
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
