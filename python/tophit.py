"""
Methods for setting up top hit workflow

get_top_hit_outputs: populates the outputs and transitions config dicts based
on the values in 'reads', 'dbs', 'filter', and 'top_alg'
"""
import os, re
from snakemake.logging import logger
from python.samples import process_sample_data
from python.qc import setup_qc_outputs
from python.annotate import get_last_alg

DEFAULT_FILTER = {'F': 0, 'B': 50}

def get_top_hit_outputs(config):
    """
    Return the list of files to be generated

    Also creates config[db_strings] that maps from database ID to file suffix
    needed to run workflow (eg: {sample}.vs.RefSeq.lastx._F0)

    Required config entries:

        dbs: dict of db names to dict containing path and format (eg: lastp)
            (currently only lastp and lastn supported)
        sample_data: dict of sample names to sample info
        filter: dict of filter_blast_m8 options (eg: {F: 0, I: 90})
        top_alg: tophit or toporg (all dbs must include tax files for taxorg)
    """
    config.setdefault('transitions', {})
    config.setdefault('outputs', set())
    db_strings = config.setdefault('db_strings', {})

    try:
        sample_data = config['sample_data']
    except KeyError:
        raise Exception("Please supply a list of samples and reads or rules "
                        "for finding them in "
                        "config[sample_data]")

    # process any patterns in sample_data[reads_patterns]
    samples = process_sample_data(sample_data, config)
    if len(samples)==0:
            raise Exception("Please supply a list of samples and reads in "
                            "config[sample_data] or rules for finding them "
                            "in config[sample_data][reads_patterns]")

    # if any samples don't have cleaned reads, try setting up QC
    if sum(1 for s in samples if 'clean' not in sample_data[s]) > 0:
        cleaned_reads_list = setup_qc_outputs(config)
        config['outputs'].update(cleaned_reads_list)
        needs_qc = True
    else:
        needs_qc = False

    try:
        dbs = config['dbs']
    except KeyError:
        raise Exception("There must be a map of databases called 'dbs' in "
                        "your configuration!")

    # setup starting point for reads for each sample
    is_prot = None
    for sample in samples:
        logger.debug('processing sample: ' + sample)
        data = sample_data[sample]
        reads_file = data['clean']
        if not isinstance(reads_file, str):
            if isinstance(reads_file[0], str):
                reads_file = reads_file[0]
            else:
                raise Exception(
                    ("I'm sorry, I only know how to work with single files "
                     " per sample. Sample {sample} seems to have multiple: "
                     "{reads_file}").format(sample=sample,
                                            reads_file=repr(reads_file)))
        extension = os.path.splitext(reads_file)[1]
        this_file_is_prot = extension == '.faa'
        if is_prot is None:
            is_prot = this_file_is_prot
        else:
            if is_prot != this_file_is_prot:
                raise Exception('All fasta files must be either nucl or AA. '
                                'some of your files have the ".faa" extension '
                                'but not all of them.')
        base_ext = extension if extension in ['.fastq', '.faa'] else '.fasta'
        logger.debug("{} becomes {}".format(extension, base_ext))
        config['transitions']['{sample}{base_ext}'.format(**vars())] = reads_file

    # set up databases
    for dbase in dbs:
        logger.debug('processing database: ' + dbase)
        search_alg = get_last_alg(dbs[dbase]['format'], is_prot)
        hit_filter = get_filter_string(config.get('filter', DEFAULT_FILTER))
        db_strings[dbase] = \
                'vs.{dbase}.{search_alg}.{hit_filter}'.format(**vars())
        topalg = config.get('top_alg', 'tophit')
        config['outputs'].add('counts.{dbase}.{topalg}.all.hitids'.format(**vars()))
        if config.get('remove_rna'):
            # also do counts without rRNA reads
            config['outputs'].add('counts.{dbase}.{topalg}.non-rRNA.hitids'.format(**vars()))
        logger.debug('added counts.{dbase}.{topalg}.hitids to outputs'.format(**vars()))

    return needs_qc

def get_filter_string(filter_dict):
    """
    return a string like _F0_B50 for a dict like {'F': 0, 'B': 50}
    """
    return "".join("_{}{}".format(k, filter_dict[k]) \
                   for k in sorted(filter_dict.keys()))

def get_non_rna_reads(wildcards, config):
    cleaned_reads = config['sample_data'][wildcards.sample]['clean']
    return re.sub(r'\.(fast[aq])$', r'.non-rRNA.\1', cleaned_reads)

def get_get_ids_cmd(wildcards, config):
    nonrna = get_non_rna_reads(wildcards, config)
    if nonrna.endswith('fasta'):
        return "grep '^>' {nonrna} | \
                 perl -pe 's/^>(\\S+).*/\\1/'".format(nonrna=nonrna)
    else:
        return "gawk '(NR+3) % 4 == 0' {nonrna} | \
                 perl -pe 's/^@(\\S+).*/\\1/'".format(nonrna=nonrna)


