"""
Methods supporting the qc workflows

setup_qc_outputs: read the sample data from config and set up the necessary
final targets (including stats files) and transitions

dummy_join_fastq: given inputs (pair of fastq files) and optional list of
sequence ids, take matching pairs from each file and combine into single
records jonied by a bunch of NNN's

"""

import re
from Bio import Seq, SeqIO, SeqRecord
from snakemake.logging import logger
from python.samples import process_sample_data
from python.tmatic import get_chemistry_barcodes

# interleaved will be automatically dropped if single file given per sample
QC_PROTOCOLS = {
    "rename": 'rename.interleaved',
    "assembly": '.'.join(['renamed',
                          'interleaved',
                          'noadapt',
                          'nophix',
                          'corrected',
                          'trimmed',
                          'dropse',
                         ]),
    "joining": '.'.join(['trim_adapt', 'joined', 'nophix']),
}

READ_DIRECTIONS = ['R1', 'R2']

def get_sample_from_reads_prefix(prefix, config):
    """
    If there is path info, assume sample is top folder name
      EG  reads/sample/reads.R1.fastq
    If there is no path info, look for _ in prefix
      EG  reads_sample.R1.fastq
    If config[sample_name] is defined, use that
    Fall back to the prefix
    """
    match = re.search(r'([^/])+/[^/]+$', prefix)
    if match:
        return match.group(1)

    match = re.search(r'^reads_(.+)$', prefix)
    if match:
        return match.group(1)

    if 'sample_name' in config:
        return config['sample_name']

    return prefix
  

def setup_qc_outputs(config):
    """
    Locate the read files for each sample and return map from sample to
    cleaned reads file.

    add 'raw' files to each sample in config[sample_data]
      (Using config[sample_data][patterns]) if missing

    Add transitions (to config[transitions]) that map
     from starting flie for QC to actual raw read files

    Add 'cleaned' file to each sample in config[sample_data]

    Works from the config dict using:
        config[sample_data][{sample}][raw]
            list of raw files in fwd, rev order
        samples_pattern:
            glob and re used to build above map (only used if map missing)
            see python.sample.collect_sample_reads
        cleaning_protocol:
            one of 'None', 'rename', 'assembly', or 'joining'

    """

    try:
        sample_data = config['sample_data']
    except KeyError:
        raise Exception("Please supply a list of samples and reads or rules "
                        "for finding them in "
                        "config[sample_data]")

    # process any patterns in sample_data[reads_patterns]
    samples = process_sample_data(sample_data)

    # find samples that need QC
    samples_with_raw_reads = [s for s in samples if 'raw' in sample_data[s]]
    if len(samples_with_raw_reads) == 0:
        raise Exception("Please supply a map from samples to raw reads "
                        " or a glob and re in"
                        " config[sample_data][reads_patterns]")

    # get protocol
    cleaning_protocol = config.get('cleaning_protocol', 'None')

    # loop back over samples and set up cleaning or interleaving if needed
    transitions = config.setdefault('transitions', {})
    outputs = []
    for sample in samples_with_raw_reads:
        raw_files = sorted(sample_data[sample]['raw'])

        # Bail out if we have too many files per sample
        if len(raw_files) > 2:
            raise Exception("I don't know how to deal with more than two"
                            " files per sample!\nSample={}\nFiles:\n{}"\
                                .format(sample,
                                        "\n".join(raw_files)))


        # starting files (define as transitions from raw files)
        if len(raw_files) == 2:
            for direction, source_file in zip(READ_DIRECTIONS, raw_files):
                transitions['{sample}.{direction}.fastq'.format(**vars())] = \
                    source_file
        else:
            transitions['{sample}.fastq'.format(**vars())] = raw_files[0]

        # cleaned suffix
        cleaned_suffix = QC_PROTOCOLS.get(cleaning_protocol, '')

        # special cases
        if cleaning_protocol != 'joining' and len(raw_files) == 1:
            # if cleaning one file for assembly, drop interleave
            cleaned_suffix = re.sub(r'\.interleaved', '', cleaned_suffix)
        if cleaning_protocol == 'None' and len(raw_files) == 2:
            # if not cleaning, but two files given, interleave them
            cleaned_suffix = 'interleaved'
        if cleaning_protocol == 'joining':
            # prepend trimmomatic params
            chemistry, barcodes = get_chemistry_barcodes(sample, config)
            cleaned_suffix = '.'.join([chemistry] + barcodes) + \
                                '.' + cleaned_suffix


        # result of QC
        cleaned_reads = '{sample}.{cleaned_suffix}.fastq'\
                                                    .format(**vars())
        sample_data[sample]['clean'] = cleaned_reads

        if config.get('remove_rna', True) in ['True', True]:
            # if rrna separation requested, add rRNA-only and non-rRNA to names
            for rrna_split in ['non-rRNA', 'rRNA-only']:
                outputs.append(re.sub(r'\.fastq$',
                                      '.{}.fastq'.format(rrna_split),
                                      cleaned_reads))
        else:
            outputs.append(cleaned_reads)

    return outputs

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
