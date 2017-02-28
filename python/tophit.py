"""
Methods for setting up top hit workflow

get_top_hit_outputs: populates the outputs and transitions config dicts based
on the values in 'reads', 'dbs', 'filter', and 'top_alg'
"""
import os
from snakemake.logging import logger

DEFAULT_FILTER = {'F': 0, 'B': 50}

def get_top_hit_outputs(config):
    """
    Return the list of files to be generated

    Also creates config[db_strings] that maps from database ID to file suffix
    needed to run workflow (eg: {sample}.vs.RefSeq.lastx._F0)

    Required config entries:

        dbs: dict of db names to dict containing path and format (eg: lastp)
            (currently only lastp and lastn supported)
        reads: dict of sample names to read files (fastq, fasta, or faa)
        filter: dict of filter_blast_m8 options (eg: {F: 0, I: 90})
        top_alg: tophit or toporg (all dbs must include tax files for taxorg)
    """
    config.setdefault('transitions', {})
    config.setdefault('outputs', set())
    db_strings = config.setdefault('db_strings', {})

    try:
        sample_reads = config['reads']
    except KeyError:
        raise Exception("There must be a map from samples to reads called "
                        "'reads' in your configuration!")

    try:
        dbs = config['dbs']
    except KeyError:
        raise Exception("There must be a map of databases called 'dbs' in "
                        "your configuration!")

    is_prot = None
    for sample, reads_file in sample_reads.items():
        logger.debug('processing sample: ' + sample)
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

    for dbase in dbs:
        logger.debug('processing database: ' + dbase)
        search_alg = get_search_alg(dbs[dbase]['format'], is_prot)
        hit_filter = get_filter_string(config.get('filter', DEFAULT_FILTER))
        db_strings[dbase] = \
                'vs.{dbase}.{search_alg}.{hit_filter}'.format(**vars())
        topalg = config.get('topalg', 'tophit')
        config['outputs'].add('counts.{dbase}.{topalg}.hitids'.format(**vars()))

def get_filter_string(filter_dict):
    """
    return a string like _F0_B50 for a dict like {'F': 0, 'B': 50}
    """
    return "".join("_{}{}".format(k, v) for k, v in filter_dict.items())

def get_search_alg(dbformat, extension):
    """
    right now looks for last db type (lastp or lastn) and extension (faa or
    not) and returns lastp, lastx, or lastn.

    Support for other dbs can be added on request.
    """
    if dbformat == 'lastp':
        if extension == 'faa':
            search_alg = 'lastp'
        else:
            search_alg = 'lastx'
    elif dbformat == 'lastn':
        if extension == 'faa':
            raise Exception("I'm sorry, I don't know how to search for faa "
                            "sequences in a lastp database!")

        else:
            search_alg = 'lastn'
    else:
        raise Exception(("I'm sorry, but the database format '{}' is not yet "
                         "supported").format(dbformat))
    return search_alg

