import re
import os
import glob
def collect_sample_reads(config, get_stats=True):
    """
    Locate the read files for each sample, set up any QC that is needed, and
    set up a 'reads' dict mapping sample names to (QCed) read files for the 
    assembly and mapping steps.

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

    reads = config.setdefault('reads',{})
    samples_pattern_data = config.setdefault('samples_pattern')
    
    # if reads alread has data and sample_pattern_data has no glob
    #  we assume that we don't want to go looking for more reads
    skip_search = len(reads) > 0 and 'glob' not in samples_pattern_data
    if not skip_search:
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

    # generate naming strings for linking QC workflow steps
    already_cleaned = samples_pattern_data.get('cleaned',True) \
                        in [True, 1, "True", "T", "true", "t"]
    if already_cleaned:
        qc_steps = []
    else:
        qc_steps = ['cleaned','corrected']
    if "min_read_length" in config:
        qc_steps.append("gte{}".format(config['min_read_length']))

    # loop back over samples and set up cleaning or interleaving if needed
    needs_qc_or_join = len(qc_steps) > 0
    transitions = config.setdefault('transitions',{})
    fasta_files = []
    for sample in list(reads.keys()):
        files = reads[sample]
        if isinstance(files, str):
            files = [files, ]
        else:
            files = sorted(files)

        # Bail out if we have too many files per sample
        if len(files)>2:
            raise Exception("I don't know how to deal with more than two"
                            " files per sample!\nSample={}\nFiles:\n{}"\
                                .format(sample,
                                        "\n".join(files)))

        # check to see if they are compressed (we can handle .gz)
        #  Bail out if one file is compressed and the other isn't
        files_gzipped = None
        for file_name in files:
            if re.search(r'\.gz$', file_name) is not None:
                if files_gzipped==False:
                    raise Exception("It seems one file is compressed and the "
                                    "other is not:\n{}".format("\n".join(files)))
                files_gzipped = True
            else:
                if files_gzipped==True:
                    raise Exception("It seems one file is compressed and the "
                                    "other is not:\n{}".format("\n".join(files)))
                files_gzipped = False

        # keep track of fasta files that will be generated
        last_fasta_file = 'reads/{sample}/reads.renamed.R12.fastq'\
                            .format(sample=sample)
        fasta_files.append(last_fasta_file)
        qc_chain = ""
        for qc_step in qc_steps:
            qc_chain += qc_step + "."
            last_fasta_file = \
                    'reads/{sample}/reads.renamed.R12.{qc_chain}fastq.gz'\
                            .format(sample=sample, qc_chain=qc_chain)
            fasta_files.append(last_fasta_file)
        
        # set target of QC, and starting point of assembly/mapping, to last file
        reads[sample] = last_fasta_file
                
        ## Create links from start of workflow to the source read files
        # are the originals gzipped?
        suffix = ".gz" if files_gzipped else ""
        # Do we have a pair of files or single file?
        if len(files)==1:
            # we are starting from interleaved, just link it in
            transitions['reads/{sample}/reads.renamed.R12.fastq{suffix}'\
                                .format(**vars())] = files[0]
        else:
            needs_qc_or_join = True
            # we are starting from paired files, link those
            transitions['reads/{sample}/reads.R1.fastq{suffix}'\
                            .format(**vars())] = files[0]
            transitions['reads/{sample}/reads.R2.fastq{suffix}'\
                            .format(**vars())] = files[1]
            
    if get_stats:
        # add stats file to final outputs for each fasta file generated
        outputs = config.setdefault('outputs', set())
        for ext in ['stats', 'hist']:
            for file_name in fasta_files:
                outputs.add(".".join((os.path.join('stats',file_name), ext)))

    return needs_qc_or_join

