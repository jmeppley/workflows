"""
Methods supporting the qc workflows

setup_qc_outputs: read the sample data from config and set up the necessary
final targets (including stats files) and transitions

dummy_join_fastq: given inputs (pair of fastq files) and optional list of
sequence ids, take matching pairs from each file and combine into single
records jonied by a bunch of NNN's

"""

import re
import glob
from collections import OrderedDict
from Bio import Seq, SeqIO, SeqRecord
from snakemake.logging import logger

def collect_sample_reads(samples_pattern_data):
    """
    Use the samples_pattern entry in config to locate read files and group into
    samples

    The samples_pattern dict should look like this:
        samples_pattern:
            glob: "../data/*.fastq"
            re: "/([^_]+)_[^/]+\.fastq"
            name: all
            metadata: samples.all.tsv

    only the 're' and 'glob' fields are used here. The glob finds files using a
    filesystem wildcard and the re should identify (as the first matched group)
    the sample name in the found file names.

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
    if len(read_files)==0:
        raise Exception(
            "The sample reads wildcard '{}' did not match any files!"\
                            .format(read_file_glob)
        )

    # collect files into lists by sample
    for read_file in read_files:
        m = sample_RE.search(read_file)
        if m is None:
            raise Exception(
                ("The sample matching expression ({}) failed to find a sample "
                 "name in the path: {}").format(sample_pattern, read_file)
            )
        sample = m.group(1)
        # sanitize sample name
        sample = re.sub(r'[^A-Za-z0-9_]','_',sample)
        reads.setdefault(sample,[]).append(read_file)

    return reads

QC_PROTOCOLS = {
    None: [],
    "assembly": ['renamed', 
                 'interleaved',
                 'noadapt',
                 'nophix',
                 'corrected',
                 'trimmed'],
    "joining": ['noadapt', 'joined', 'trimmed', 'nophix'],
}

STEP_FASTQ_OUTPUTS = {
    'noadapt'

QC_STEP_MERGES_PAIRS = set('interleaved', 'joined')

READ_DIRECTIONS = ['R1', 'R2']

def get_qc_steps(config):
    """
    return list of file suffixes necessary to implement selected protocol
    """
    cleaning_protocol = config.get('cleaning_protocol', None)
    qc_steps = QC_PROTOCOLS[cleaning_protocol]
    return qc_steps

def build_qc_chain(sample, paired, steps,
                   fasta_files,
                   cleand_reads_by_sample):
    """
    given the sample and the steps chain
    """
    # keep track of fasta files that will be generated
    #  starting with first file in chain and gothrough qc_steps
    if paired:
        starting_files = [
            'reads/{sample}/reads.{dir}.fastq'.format(**vars()) \
                for dir in READ_DIRECTIONS
        ]
    else:
        starting_files = [
            'reads/{sample}/reads.fastq'.format(**vars()),
        ]
    fasta_files.extend(starting_files)

    # want to keep track of all fasta files, start with starting files
    last_fasta_files = starting_files
    for step in steps:
        last_fasta_files = get_next_fasta_files(last_fasta_files,
                                                 step,
                                                 paired)
        # no longer paired if we went through a merge step
        if len(last_fasta_files)==1:
            paired = False

    # we should end with one file
    if paired:
        raise Exception((
            "end of QC chain must end in single (interleaved or "
            "merged) file. Sample {} has 2 files, but there is no "
            "merging step in chain: {}".format(
                sample,
                repr(steps)))
    cleaned_reads[sample] = last_fasta_files[0]

    return starting_files


def get_next_fasta_files(last_fasta_files,
                         step,
                         paired):
    """
    Given a starting file(s) and a qc step, get next file or files
    """
    if step not in QC_STEP_MERGES_PAIRS:
        if paired:
            old_last_fasta_files = list(last_fasta_files)
            last_fasta_files = [
                new_last_fasta_files.append(
                    get_next_fasta_file(last_fasta_file,
                                        step,
                                        fasta_files)
                ) \
                for last_fasta_file in last_fasta_files
            ]
        else:
            last_fasta_files = [
                new_last_fasta_file(last_fasta_files[0],
                                    step,
                                    fasta_files)
            ]
    else:
        if not paired:
            raise Exception((
                "The qc protocol {} with steps {} includes a merging step "
                "({}), but sample {} only has one file. We cannot "
                "make sense of this. Sorry") \
              .format(config.get('cleaning_protocol',None),
                       repr(steps),
                       step,
                       sample)
            )
        last_fasta_files = [
            new_last_fasta_file_merged(last_fasta_files,
                                       step,
                                       fasta_files)
        ]
    return last_fasta_files


def new_last_fasta_file_merged(last_fasta_files, step, fasta_files):
    """
    drops the R1 from the file name before processing
    """
    dropped_dir = re.sub(r'\.R1', '', last_fasta_files[0])
    return new_last_fasta_file(dropped_dir, step, fasta_files)


def new_last_fasta_file(last_fasta_file, step, fasta_files):
    """
    Given starting file(s) and step name
    
     * figure out the next file step in chain 
     * add all intermediate fasta files to list
    """

    # by default, insert step name before .fastq or .R1.fastq
    default_rexp = r'(?<!\.R[12])((?:\.R1)?\.fastq)'
    defualt_subst = r'.{step}\1'.format(step=step)
    default_fastq_outputs = (default_rexp, defualt_subst)

    # get step specific substitutions
    fastq_outputs = STEP_FASTQ_OUTPUTS.get(step,
                                           default_fastq_outputs)

    # use substitutions to populate fasta file list
    for rexp, subst in fastq_outputs:
        fasta_files.append(re.sub(rexp, subst, 
                                  last_fasta_file))

    # use default subst to build up final file name
    #  (this particular file name may never get made)
    next_fasta_file = re.sub(default_rexp,
                             default_subst,
                             last_fasta_file)
    return next_fasta_file

def setup_qc_outputs(config, get_stats=True):
    """
    Locate the read files for each sample, set up any QC that is needed, and
    set up a 'cleaned_reads' dict mapping sample names to (QCed) read files 
    for the assembly and mapping steps.

    If there is not already a "reads" dict in config, build it from the samples
    pattern. (see collect_sample_reads() above)

    returns list of snakefiles to include to get necessary QC rules

    QC rules depend on config['cleaning_protocol'] which can be one of:
        assembly (alias for assembly-bfc)
        assembly-bfc 
        joining (alias for joing-pandaseq)
        joining-pandaseq
        joining-pear
        joining-flash

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
    
    # if reads already has data and sample_pattern_data has no glob
    #  we assume that we don't want to go looking for more reads
    #  (we have defaults to try if nothing supplied)
    if len(reads) == 0 or 'glob' in samples_pattern_data:
        reads.update(collect_sample_reads(samples_pattern_data))

    # generate naming strings for linking QC workflow steps
    qc_steps = get_qc_steps(config)

    # loop back over samples and set up cleaning or interleaving if needed
    snakefiles = []
    transitions = config.setdefault('transitions',{})
    fasta_files = []
    cleaned_reads = {}
    for sample in list(reads.keys()):
        files = sorted(reads[sample])

        # Bail out if we have too many files per sample
        if len(files)>2:
            raise Exception("I don't know how to deal with more than two"
                            " files per sample!\nSample={}\nFiles:\n{}"\
                                .format(sample,
                                        "\n".join(files)))

        starting_files = build_qc_chain(sample, qc_steps, len(files) == 2,
                                        fasta_files,
                                        cleaned_reads)


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
                outputs.add(".".join((file_name, ext)))

    return needs_qc_or_join


def rev_comp_rec(record, qual=False, suffix=''):
    """
    reverse complement the sequence of a SeqRecord object
    """
    new_record = SeqRecord.SeqRecord(record.seq.reverse_complement(),
                                     id=record.id+suffix,
                                     name=record.name,
                                     description=record.description)
    if qual:
        new_record.letter_annotations['phred_quality'] = \
                list(reversed(record.letter_annotations['phred_quality']))
    return new_record

READ_INDEX_DIR_RE = re.compile(r'(#(\d+|[ACTGN]+))?(/[12])?')
DUP_COLON_RE = re.compile(r':+')
def sanitize_record_id_set(record_id_set):
    """ make sure collection is a set and remove runs of colons """
    prev_len = len(record_id_set)
    record_id_set = set(READ_INDEX_DIR_RE.sub('', DUP_COLON_RE.sub(r':', s)) \
                            for s in record_id_set)

    if prev_len != len(record_id_set):
        ## Uh oh, this hack broke something
        raise Exception("The list of unpaired reads has name collisions "
                        "if multiple colons are ignored. Your reads are "
                        "using a naming convetion that the author of this "
                        "software (illuminPrep) didn't anticipate. My "
                        "apologies")

    return record_id_set

def sanitize_record_id(record):
    """
    remove direction and index suffix and any strings of colons
    """
    record.id = READ_INDEX_DIR_RE.sub('', DUP_COLON_RE.sub(r':', record.id))

class EndOfRecords(Exception):
    """
    Special exception to indicate that we reached the end of a record iterator
    while processing it inside a generator.
    """

def get_next_pair(rec_iter, rec_cache, pair_cache, is_rev=False):
    """
    Get the next record from iterator
    If it has no pair in the pair cache, add it to rec_cache and quit
    If it does have a pair
        yield every rec in cache with null pair
        yield rec and pair
    """
    try:
        record = next(rec_iter)
    except StopIteration:
        # special exception to differenticate from end of generator
        raise EndOfRecords()
    else:
        sanitize_record_id(record)
        try:
            pair = pair_cache.pop(record.id)
        except KeyError:
            # We have not encountered the rev, yet, save it
            rec_cache[record.id] = record
        else:
            # We found a matching record in pair cache
            while len(rec_cache) > 0:
                # assume files in order, so prev records have no pair
                if is_rev:
                    yield (None, rec_cache.popitem()[1])
                else:
                    yield (rec_cache.popitem()[1], None)
            # now yield the pair we just found
            yield record, pair

def merge_record_iters(fwd_iter, rev_iter, record_id_filter=None):
    """
    Iterate over a pair of sequences files with the same order of read names
    but allow for missing records in each
    """
    if record_id_filter is not None:
        record_id_filter = sanitize_record_id_set(record_id_filter)

    fwd_record_cache = {}
    rev_record_cache = {}
    reached_both_ends = False
    while not reached_both_ends:

        # get next fwd read and see if we've already seen its pair
        try:
            for pair in get_next_pair(fwd_iter,
                                      fwd_record_cache,
                                      rev_record_cache):
                yield pair
        except EndOfRecords:
            end_of_fwd_records = True
        else:
            end_of_fwd_records = False

        # get next rev read and see if we've already seen its pair
        try:
            for pair in get_next_pair(rev_iter,
                                      rev_record_cache,
                                      fwd_record_cache,
                                      is_rev=True):
                yield pair
        except EndOfRecords:
            end_of_rev_records = True
        else:
            end_of_rev_records = False

        reached_both_ends = end_of_fwd_records and end_of_rev_records

    # empty caches
    while len(fwd_record_cache) > 0:
        yield (fwd_record_cache.popitem()[1], None)
    while len(rev_record_cache) > 0:
        yield (None, rev_record_cache.popitem()[1])

def dummy_join_fastq(inputs,
                     outputs,
                     log_files,
                     batch_size=10000,
                     gap=20,
                     **kwargs
                    ):
    """
    Join reads from fwd and rev fastq files with a string of {gap} N's

     * inputs is class from snake rule with "fwd" and "rev" file names
     * outputs is list from snake rule with output file as first element
     * log_files is list with log file as first element

    Paramters:
        
     * If record_id_filter given, only join read pairs listed there.
     * batch_size is how many output records to cache between writes
     * gap is how many Ns' to insert between ends

    """
    if gap >= 0:
        gap_seq = 'N'*gap
        gap_qual = [0]*gap

    counts = {
        'join_count': 0,
        'fwd_count': 0,
        'rev_count': 0,
        'total_joined': 0,
        'total_written': 0,
    }
    faked_joins = []

    fwd_records = SeqIO.parse(inputs.fwd, 'fastq')
    rev_records = SeqIO.parse(inputs.rev, 'fastq')
    with open(outputs[0], 'w') as out_fastq_stream:
        for frec, rrec in merge_record_iters(fwd_records, rev_records,
                                             **kwargs):
            # join seqs
            new_records = []
            if frec is None and rrec is None:
                logger.warning("Both ends missing from input") # this shouldn't
                continue
            if frec is None:
                logger.debug("Forward seq trimmed to oblivion")
                new_records.append(rev_comp_rec(rrec, qual=True, suffix=".rev"))
                counts['rev_count'] += 1
            elif rrec is None:
                logger.debug("Reverse seq trimmed to oblivion")
                new_records.append(frec)
                counts['fwd_count'] += 1
            elif gap >= 0:
                counts['join_count'] += 1
                # join sequence
                new_seq = frec.seq \
                          + Seq.Seq(gap_seq, frec.seq.alphabet) \
                          + rrec.seq.reverse_complement()
                new_record = SeqRecord.SeqRecord(new_seq,
                                                 id=frec.id,
                                                 name=frec.name,
                                                 description="Faked join")
                # join quality
                new_record.letter_annotations['phred_quality'] = \
                        frec.letter_annotations['phred_quality'] + \
                        gap_qual + \
                        list(reversed(rrec.letter_annotations['phred_quality']))
                new_records.append(new_record)
            else:
                # gap < 0 means don't join...add separately
                new_records.append(frec)
                new_records.append(rev_comp_rec(rrec, qual=True, suffix=".rev"))

            faked_joins.extend(new_records)
            if len(faked_joins) >= batch_size:
                n_written = SeqIO.write(faked_joins,
                                        out_fastq_stream,
                                        format='fastq')
                if n_written != len(faked_joins):
                    logger.warning("Only %d of %d faked joins written!" %
                                   (n_written, len(faked_joins)))
                counts['total_joined'] += len(faked_joins)
                counts['total_written'] += n_written
                del faked_joins[:]

        # at end of loop, write remaining cached records
        n_written = SeqIO.write(faked_joins, out_fastq_stream, format='fastq')
        if n_written != len(faked_joins):
            logger.warning("Only %d of %d faked joins written!" %
                           (n_written, len(faked_joins)))
        counts['total_joined'] += len(faked_joins)
        counts['total_written'] += n_written

        # Report some counts
        msg = """
#======================
# Faked joins
#  Total written: {total_written} of {total_joined}
#  Dummy Joins: {join_count}
#  FwdOnly: {fwd_count}
#  RevOnly: {rev_count}
#======================
""".format(**counts)

        with open(log_files[0], 'a') as log_out_stream:
            log_out_stream.write(msg)

        logger.debug(msg)
