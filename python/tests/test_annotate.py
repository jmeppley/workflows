"""
Tests for annotate.py
"""
import re
import yaml
from python.annotate import \
        get_batch_launcher_prefix, \
        get_db_dot_fmt_strings, \
        get_db_types, \
        get_last_alg

class Wildcards():
    """ Dummy class to take the place of wildcards in tests """
    def __init__(self,):
        pass

def test_batch_launcher_prefix():
    """ Should return empty string if no frags """
    wildcards = Wildcards()
    wildcards.db = 'test_db_1'
    config = {'dbs': {'test_db_1': {'frags': 1},
                      'test_db_2': {'frags': 4},
                     }}
    prefix = get_batch_launcher_prefix(wildcards, config)
    assert prefix == ''

    wildcards = Wildcards()
    wildcards.db = 'test_db_2'
    prefix = get_batch_launcher_prefix(wildcards, config)
    assert re.search(r' -N 4 ', prefix) is not None


def test_get_last_alg():
    """ make search suffix reflect  query type """
    assert get_last_alg('lastp', 'faa') == 'lastp'
    assert get_last_alg('lastp', 'ffn') == 'lastx'
    assert get_last_alg('lastp', 'fasta') == 'lastx'
    assert get_last_alg('lastn', 'fna') == 'lastn'
    assert get_last_alg('bwadb', 'fasta') == 'bwa.sam'

def test_dbs():
    """
    Use the two included sample configs as test cases for
    the get_db_types, and get_db_dot_fnt_strings methods
    """
    with open('annotation.yaml.sample') as config_handle:
        config = yaml.load(config_handle)
    gene_family_dbs, tax_db = get_db_types(config)
    assert gene_family_dbs == ['COG',]
    assert tax_db == 'RefSeq'
    db_fmt_strings = get_db_dot_fmt_strings(gene_family_dbs,
                                            config,)
    assert db_fmt_strings == ['COG.tbl']

    with open('docker/config.yaml') as config_handle:
        config = yaml.load(config_handle)
    gene_family_dbs, tax_db = get_db_types(config)
    assert sorted(gene_family_dbs) == ['COG', 'KEGG', 'PFAM', 'TIGR']
    assert tax_db == 'RefSeq'
    db_fmt_strings = get_db_dot_fmt_strings(gene_family_dbs,
                                            config,
                                            'fna')
    assert sorted(db_fmt_strings) == ['COG.tbl',
                                      'KEGG.lastx',
                                      'PFAM.tbl',
                                      'TIGR.tbl',
                                     ]

