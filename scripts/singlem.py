import logging
import os
import re
from python.download import remote_wrapper

def get_sample_data(config):
    """
    Look at config['sample_data'] (should be a dict) and
    wildcard glob strings in config['reads_glob'] config['contigs_glob']
    
    sample list is the union of explicit samples in sample_data and 
    anything matched by reads_glob. Filse found with contigs_glob
    are only used if sample aready selected. 

    Sample list from over-inclusive reads_globs can be pared down with
    config['sample_filter_rexp']
    
    """

    # optionally filter samples
    if 'sample_filter_rexp' in config:
        sf_rexp = re.compile(config['sample_filter_rexp'])
        def sample_filter(sample):
            return sf_rexp.search(sample) is not None
    else:
        def sample_filter(sample):
            return True

    # check if any hard coded files need to be downloaded
    sample_data = {s:{t:remote_wrapper(f, config) \
                      for t,f in files.items()} \
                   for s,files in config.get('sample_data', {}).items()}

    # look for reads with wildcard
    if "reads_glob" in config:
        glob_sample_data = {}
        # string like "/path/to/files/{sample}.fastq"
        reads_glob = config['reads_glob']

        wildcard_values = remote_wrapper(reads_glob, config, glob=True)
        
        # which wildcard in glob was "sample"
        wildcard_names = list(wildcard_values._fields)
        sample_index = wildcard_names.index('sample')
        
        # build sample_data dictionary
        for wildcard_tuple in zip(*list(wildcard_values)):
            fmt_dict = dict(zip(wildcard_names, wildcard_tuple))
            logging.debug(fmt_dict)
            for value in wildcard_tuple:
                # skip if wildcard spans folders
                if re.search(os.path.sep, value):
                    break
            else:
                # add sample
                sample = wildcard_tuple[sample_index]
                if not sample_filter(sample):
                    continue
                glob_sample_data[sample] = \
                        {'reads':
                         remote_wrapper(reads_glob.format(**fmt_dict), config)}

        if len(glob_sample_data) == 0:
            logger.warning("Wildcard glob DID NOT match any files: " +
                           reads_glob)

        else:
            sample_data.update(glob_sample_data)

    if len(sample_data) == 0:
        raise Exception("Please supply either a map of sample names to "
                        "read files in 'sample_data' or a wildcard glob "
                        "string in 'reads_glob' to locate files. "
                        "EG: reads_glob=/path/to/files/{sample}.fastq")

    # look for contigs with wildcard
    if "contigs_glob" in config:
        # look for all matches and filter, so we can use remote provider
        contigs_glob = config['contigs_glob']
        wildcard_values = remote_wrapper(contigs_glob, config, glob=True)
        
        # which wildcard in glob was "sample"
        wildcard_names = list(wildcard_values._fields)
        sample_index = wildcard_names.index('sample')
 
        # update sample_data dictionary
        for wildcard_tuple in zip(*list(wildcard_values)):
            fmt_dict = dict(zip(wildcard_names, wildcard_tuple))
            logging.debug(fmt_dict)
            for value in wildcard_tuple:
                # skip if wildcard spans folders
                if re.search(os.path.sep, value):
                    break
            else:
                # add sample
                sample = wildcard_tuple[sample_index]
                if sample not in sample_data:
                    # we're only interested in things we have reads for
                    continue
                if 'contigs' in sample_data[sample]:
                    # don't overwrite hard coded files
                    continue
                sample_data[sample]['contigs'] = \
                     remote_wrapper(contigs_glob.format(**fmt_dict), config)

    config['sample_data'] = sample_data
    return sample_data

