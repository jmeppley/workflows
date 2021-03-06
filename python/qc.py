"""
Methods supporting the qc workflows

setup_qc_outputs: read the sample data from config and set up the necessary
final targets (including stats files) and transitions

dummy_join_fastq: given inputs (pair of fastq files) and optional list of
sequence ids, take matching pairs from each file and combine into single
records jonied by a bunch of NNN's

"""

import re
from collections import defaultdict
import yaml
from Bio import Seq, SeqIO, SeqRecord
from snakemake.logging import logger
from python.samples import process_sample_data
from python.tmatic import get_chemistry_barcodes
from python.common import get_file_name

# from python.common import is_in_working_dir

# interleaved will be automatically dropped if single file given per sample
QC_PROTOCOLS = {
    "rename": "rename.interleaved",
    "interleave": "interleaved",
    "ec_only": "corrected.dropse",
    "assembly_fastp": ".".join(
        ["renamed", "interleaved", "cleaned", "corrected", "dropse",]
    ),
    "assembly_fastp_no_ec": ".".join(["renamed", "interleaved", "cleaned", "dropse",]),
    "assembly": ".".join(
        [
            "renamed",
            "interleaved",
            "noadapt",
            "nophix",
            "corrected",
            "trimmed",
            "dropse",
        ]
    ),
    "anvio": "renamed-QUALITY_PASSED_interleaved",
    "anvio_pe": "renamed-QUALITY_PASSED_interleaved",
    "assembly_no_ec": ".".join(
        ["renamed", "interleaved", "noadapt", "nophix", "trimmed", "dropse",]
    ),
    "joining": ".".join(["trim_adapt", "joined", "nophix"]),
    "join_and_standards": ".".join(["trim_adapt", "joined", "nophix", "nospikes"]),
}

# list of protocols with no error correction
NON_EC_PROTOCOLS = [
    "joining",
    "assembly_no_ec",
    "assembly_fastp_no_ec",
    "join_and_standards",
]
JOINING_PROTOCOLS = ["joining", "join_and_standards"]

READ_DIRECTIONS = ["R1", "R2"]


def get_sample_from_reads_prefix(prefix, config):
    """
    check for two specific cases, then fall back to the prefix:
    If there is path info, assume sample is top folder name
      EG  reads/sample/reads.R1.fastq
    If there is no path info, look for reads_ in prefix
      EG  reads_sample.R1.fastq
    Fall back to the prefix
    """
    match = re.search(r"([^/])+/reads$", prefix)
    if match:
        return match.group(1)

    match = re.search(r"^reads_(.+)$", prefix)
    if match:
        return match.group(1)

    return prefix


def setup_qc_outputs(config):
    """
    Locate the read files for each sample and return map from sample to
    cleaned reads file.

    add 'raw' files to each sample in config[sample_data]
      (Using config[sample_data][reads_patterns]) if missing

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
        sample_data = config["sample_data"]
    except KeyError:
        raise Exception(
            "Please supply a list of samples and reads or rules "
            "for finding them in "
            "config[sample_data]"
        )

    # process any patterns in sample_data[reads_patterns]
    samples = process_sample_data(sample_data, config)

    # initialize:
    #  create symlinks between paired files
    transitions = config.setdefault("transitions", {})
    #  dict of filter transitions
    filter_dict = config.setdefault("filter_dict", {})
    #  list of cleaned reads to include in assembly
    cleaned_read_list = config.setdefault("cleaned_read_list", [])
    #  list of final read files to consider as outputs
    outputs = []

    # LOOP over samples
    for sample in samples:
        data = sample_data[sample]

        ## CLEAN
        # Did this sample come with clean reads? If so:
        #  1) set up pairs to be interleaved
        #  2) replace with locally named file and add transition if not in workdir
        if "clean" in data:
            cleaned_reads = process_cleaned_reads(sample, data, transitions)
            if cleaned_reads != None:
                cleaned_read_list.append(cleaned_reads)

        ## RAW
        # Does this sample have raw reads that need QC?
        if "raw" in data and "clean" not in data:

            # get files as sorted list
            raw_files = (
                [data["raw"],] if isinstance(data["raw"], str) else sorted(data["raw"])
            )

            # get protocol (try sample_Data first, fall back to global)
            cleaning_protocol = data.setdefault(
                "protocol", config.setdefault("cleaning_protocol", "None")
            )

            # check to see if input is compressed (we can handle .gz)
            extension = "fastq.gz" if are_files_gzipped(raw_files) else "fastq"

            # do we need to filter reads first?
            #  WARNING: filtering only works with uncompressed fastq
            filter_file = data.get("filter", None)

            # starting files (define as transitions from raw files)
            new_raw_files = []
            if len(raw_files) > 2:
                # The reads are probably split up into lanes, merge into one pair
                for direction, file_list in setup_merge_by_lanes(
                    raw_files, sample
                ).items():
                    local_raw_file = "{sample}.{direction}.{extension}".format(**vars())
                    filter_or_link_transition(
                        file_list, local_raw_file, filter_file, transitions, filter_dict
                    )
                    new_raw_files.append(local_raw_file)
            elif len(raw_files) == 2:
                for direction, source_file in zip(READ_DIRECTIONS, raw_files):
                    local_raw_file = "{sample}.{direction}.{extension}".format(**vars())
                    filter_or_link_transition(
                        source_file,
                        local_raw_file,
                        filter_file,
                        transitions,
                        filter_dict,
                    )
                    new_raw_files.append(local_raw_file)
            else:
                local_raw_file = "{sample}.{extension}".format(**vars())
                filter_or_link_transition(
                    raw_files[0], local_raw_file, filter_file, transitions, filter_dict
                )
                new_raw_files.append(local_raw_file)

            # everything else should work from the linked files
            data["raw"] = new_raw_files

            # cleaned suffix
            cleaned_suffix = get_cleaned_suffix(cleaning_protocol, sample, config)

            ## result of QC
            # the actual QC output depends on the remove_rna flag
            if config.get("remove_rna", True) in ["True", True]:
                # if rrna separation requested . . .
                logger.debug("adding separated rrna reads to output")
                # add rRNA reads to output
                for rrna_split in ["SSU", "LSU"]:
                    outputs.append(
                        "{sample}.{cleaned_suffix}{rrna_split}.fastq".format(**vars())
                    )
                # set non-rRNA as cleaned target
                cleaned_raw_reads = "{sample}.{cleaned_suffix}non-rRNA.fastq".format(
                    **vars()
                )
            else:
                # otherwise, QC ends at cleaned_suffix
                cleaned_raw_reads = "{sample}.{cleaned_suffix}fastq".format(**vars())

            # this is the alias that will point to the actual QC output
            cleaned_reads = "{sample}.clean.fastq".format(**vars())
            cleaned_read_list.append(cleaned_reads)
            outputs.append(cleaned_reads)
            transitions[cleaned_reads] = cleaned_raw_reads
            data["clean"] = cleaned_reads

            # samples.get_sample_reads_for_mapping needs this
            data["clean_realpath"] = cleaned_raw_reads

    return outputs


def filter_or_link_transition(
    source_file, dest_file, filter_file, transitions, filter_dict,
):
    """ Both the transitions mechanism and the filter mechanism can
    handle multiple files """
    if filter_file is not None:
        filter_dict[dest_file] = {"all_reads": source_file, "filter": filter_file}
    else:
        transitions[dest_file] = source_file


def process_cleaned_reads(sample, data, transitions):
    """
    For a sample that came with clean reads, make sure we end up with
    a single interleaved file of reads in working folder
    """
    drop_clean_reads = False
    remote_cleaned_reads = data["clean"]

    # send reads through QC if there is any processing that neds to be done:
    #  EG: interleaving or filtering
    needs_filter = "filter" in data
    needs_interleaving = (
        not isinstance(remote_cleaned_reads, str) and len(remote_cleaned_reads) > 1
    )
    if needs_filter or needs_interleaving:
        # we do this by sending through QC with no protocol
        data.setdefault("raw", remote_cleaned_reads)
        data.setdefault("protocol", "None")
        # remove 'clean' so QC will get setup below
        del data["clean"]
        return None
    elif len(remote_cleaned_reads) == 0:
        raise Exception("No clean reads for sample {}:\n{}".format(sample, repr(data)))
    else:
        # change list with one string into single str object
        if not isinstance(remote_cleaned_reads, str):
            remote_cleaned_reads = remote_cleaned_reads[0]
        # build standardized cleaned reads file name in workdir
        local_cleaned_reads = "{sample}.clean.fastq".format(**vars())
        if re.search(r"\.gz$", remote_cleaned_reads) is not None:
            local_cleaned_reads += ".gz"
        # if input file does not conform to naming convention, create link
        if local_cleaned_reads != remote_cleaned_reads:
            transitions[local_cleaned_reads] = remote_cleaned_reads
        data["clean"] = local_cleaned_reads
        return local_cleaned_reads


def are_files_gzipped(raw_files):
    """ if all input files end in .gz, return true
        throw exception if mismatched """
    files_are_gzipped = None
    for file_name in raw_files:
        if re.search(r"\.gz$", file_name) is not None:
            if files_are_gzipped is False:
                raise Exception(
                    "It seems one file is compressed and the "
                    "other is "
                    "not:\n{}".format("\n".join(raw_files))
                )
            files_are_gzipped = True
        else:
            if files_are_gzipped:
                raise Exception(
                    "It seems one file is compressed and the "
                    "other is "
                    "not:\n{}".format("\n".join(raw_files))
                )
            files_are_gzipped = False
    return files_are_gzipped


def get_cleaned_suffix(cleaning_protocol, sample, config):
    """ determine what the final output file name suffix is """
    cleaned_suffix = QC_PROTOCOLS.get(cleaning_protocol, "")
    num_raw_files = len(config["sample_data"][sample]["raw"])

    # special cases
    if cleaning_protocol not in JOINING_PROTOCOLS and num_raw_files == 1:
        # if cleaning one file for assembly, drop interleave
        cleaned_suffix = re.sub(r"interleaved", "", cleaned_suffix)
    if cleaning_protocol == "None" and num_raw_files == 2:
        # if not cleaning, but two files given, interleave them
        cleaned_suffix = "interleaved"
    if cleaning_protocol in JOINING_PROTOCOLS:
        # prepend trimmomatic params
        chemistry, barcodes = get_chemistry_barcodes(sample, config)
        barcodes = ".".join(barcodes)
        cleaned_suffix = ".".join([chemistry, barcodes]) + "." + cleaned_suffix

    # make sure it ends in a dot
    cleaned_suffix += "."

    # remove dupmlicate dots left over from above tweaks
    cleaned_suffix = re.sub(r"\.\.", ".", cleaned_suffix)

    return cleaned_suffix


def setup_merge_by_lanes(raw_files, sample):
    """
    Pair files, and setup merge into a single pair of files
    """
    pairs = defaultdict(dict)
    for file_name in raw_files:
        # key off of file name with "R1" or "R2" removed
        if re.search(r"[-._]R1[-._]", file_name):
            pair_key = re.sub(r"[-._]R1[-._]", "", file_name)
            pairs[pair_key]["R1"] = file_name
        elif re.search(r"[-._]R2[-._]", file_name):
            pair_key = re.sub(r"[-._]R2[-._]", "", file_name)
            pairs[pair_key]["R2"] = file_name
        else:
            raise Exception("This does not appear to be a pairred file! " + file_name)

    logger.debug(yaml.dump(pairs))

    concatenations = {"R1": [], "R2": []}
    for pair in sorted(pairs.keys()):
        for direction in concatenations.keys():
            if direction not in pairs[pair]:
                raise Exception("Lane has no pair: " + repr(pairs[pair]))
            concatenations[direction].append(pairs[pair][direction])

    logger.debug(yaml.dump(concatenations))
    return concatenations


def rev_comp_rec(record, qual=False, suffix=""):
    """
    reverse complement the sequence of a SeqRecord object
    """
    new_record = SeqRecord.SeqRecord(
        record.seq.reverse_complement(),
        id=record.id + suffix,
        name=record.name,
        description=record.description,
    )
    if qual:
        new_record.letter_annotations["phred_quality"] = list(
            reversed(record.letter_annotations["phred_quality"])
        )
    return new_record


READ_INDEX_DIR_RE = re.compile(r"(#(\d+|[ACTGN]+))?(/[12])?")
DUP_COLON_RE = re.compile(r":+")


def sanitize_record_id_set(record_id_set):
    """ make sure collection is a set and remove runs of colons """
    prev_len = len(record_id_set)
    record_id_set = set(
        READ_INDEX_DIR_RE.sub("", DUP_COLON_RE.sub(r":", s)) for s in record_id_set
    )

    if prev_len != len(record_id_set):
        ## Uh oh, this hack broke something
        raise Exception(
            "The list of unpaired reads has name collisions "
            "if multiple colons are ignored. Your reads are "
            "using a naming convetion that the author of this "
            "software (illuminPrep) didn't anticipate. My "
            "apologies"
        )

    return record_id_set


def sanitize_record_id(record):
    """
    remove direction and index suffix and any strings of colons
    """
    record.id = READ_INDEX_DIR_RE.sub("", DUP_COLON_RE.sub(r":", record.id))


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
            for pair in get_next_pair(fwd_iter, fwd_record_cache, rev_record_cache):
                yield pair
        except EndOfRecords:
            end_of_fwd_records = True
        else:
            end_of_fwd_records = False

        # get next rev read and see if we've already seen its pair
        try:
            for pair in get_next_pair(
                rev_iter, rev_record_cache, fwd_record_cache, is_rev=True
            ):
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


def dummy_join_fastq(inputs, outputs, log_files, batch_size=10000, gap=20, **kwargs):
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
        gap_seq = "N" * gap
        gap_qual = [0] * gap

    counts = {
        "join_count": 0,
        "fwd_count": 0,
        "rev_count": 0,
        "total_joined": 0,
        "total_written": 0,
    }
    faked_joins = []

    fwd_records = SeqIO.parse(get_file_name(inputs.fwd), "fastq")
    rev_records = SeqIO.parse(get_file_name(inputs.rev), "fastq")
    with open(get_file_name(outputs), "w") as out_fastq_stream:
        for frec, rrec in merge_record_iters(fwd_records, rev_records, **kwargs):
            # join seqs
            new_records = []
            if frec is None and rrec is None:
                logger.warning("Both ends missing from input")  # this shouldn't
                continue
            if frec is None:
                logger.debug("Forward seq trimmed to oblivion")
                new_records.append(rev_comp_rec(rrec, qual=True, suffix=".rev"))
                counts["rev_count"] += 1
            elif rrec is None:
                logger.debug("Reverse seq trimmed to oblivion")
                new_records.append(frec)
                counts["fwd_count"] += 1
            elif gap >= 0:
                counts["join_count"] += 1
                # join sequence
                new_seq = (
                    frec.seq
                    + Seq.Seq(gap_seq, frec.seq.alphabet)
                    + rrec.seq.reverse_complement()
                )
                new_record = SeqRecord.SeqRecord(
                    new_seq, id=frec.id, name=frec.name, description="Faked join"
                )
                # join quality
                new_record.letter_annotations["phred_quality"] = (
                    frec.letter_annotations["phred_quality"]
                    + gap_qual
                    + list(reversed(rrec.letter_annotations["phred_quality"]))
                )
                new_records.append(new_record)
            else:
                # gap < 0 means don't join...add separately
                new_records.append(frec)
                new_records.append(rev_comp_rec(rrec, qual=True, suffix=".rev"))

            faked_joins.extend(new_records)
            if len(faked_joins) >= batch_size:
                n_written = SeqIO.write(faked_joins, out_fastq_stream, format="fastq")
                if n_written != len(faked_joins):
                    logger.warning(
                        "Only %d of %d faked joins written!"
                        % (n_written, len(faked_joins))
                    )
                counts["total_joined"] += len(faked_joins)
                counts["total_written"] += n_written
                del faked_joins[:]

        # at end of loop, write remaining cached records
        n_written = SeqIO.write(faked_joins, out_fastq_stream, format="fastq")
        if n_written != len(faked_joins):
            logger.warning(
                "Only %d of %d faked joins written!" % (n_written, len(faked_joins))
            )
        counts["total_joined"] += len(faked_joins)
        counts["total_written"] += n_written

        # Report some counts
        msg = """
#======================
# Faked joins
#  Total written: {total_written} of {total_joined}
#  Dummy Joins: {join_count}
#  FwdOnly: {fwd_count}
#  RevOnly: {rev_count}
#======================
""".format(
            **counts
        )

        with open(log_files[0], "a") as log_out_stream:
            log_out_stream.write(msg)

        logger.debug(msg)
