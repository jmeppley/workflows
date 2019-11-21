"""
Functions used in the read_mapping workflows

check_sample_data(config):
    A simpler approach to getting sample data from wildcards in the
    configuration than used  by most of the other workflows. This just
    expands a wildcard glob.
"""
import re
import os
from jme.dynamic_remote_snake.remote import remote_wrapper

def check_sample_data(config):
    """ if sample_data not explicitly listed in config, work from glob """
    if "sample_glob" in config:
        # string like "/path/to/files/{sample}.fastq"
        sample_glob = config['sample_glob']
        wildcard_values = remote_wrapper(sample_glob, config, glob=True)
        sample_data = config.setdefault('sample_data', {})

        # which wildcard in glob was "sample"
        wildcard_names = list(wildcard_values._fields)
        sample_index = wildcard_names.index('sample')

        # build sample_data dictionary
        for wildcard_tuple in zip(*list(wildcard_values)):
            for value in wildcard_tuple:
                # skip if wildcard spans folders
                if re.search(os.path.sep, value):
                    break
            else:
                # add sample
                sample = wildcard_tuple[sample_index]
                sample_data[sample] = sample_glob.format(
                                            **dict(zip(wildcard_names,
                                                       wildcard_tuple)))
    if "sample_data" not in config:
            raise Exception("Please supply read files explicitly in "
                            "config[sample_data] or with  "
                            "config[sample_glob]")

    if len(config["sample_data"]) == 0:
            raise Exception("wildcard glob matched 0 files!")


    # apply remote_wrapper to all files in sample_data
    config['sample_data'] = {s:remote_wrapper(f, config) \
                             for s, f in config['sample_data'].items()}
