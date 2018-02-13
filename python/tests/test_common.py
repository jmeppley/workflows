"""
Tests for common.py
"""
import re
from python.common import get_version, apply_defaults

def test_get_version():
    """ Try a few programs """

    # cmsearch should have "INFERNAL" in the version string
    result = get_version('cmsearch', '-h', lines=1)
    match = re.search(r'INFERNAL', result)
    assert match is not None

    # samtools
    result = get_version('samtools', lines=[0,])
    match = re.search(r'samtools', result)
    assert match is not None

    # centrifuge
    result = get_version('centrifuge', lines=[0,],
                         regular_expression=r'([Vv]ersion.+)$')
    match = re.search(r'\d+\.\d+\.\d+', result)
    assert match is not None

    # lastal
    result = get_version('lastal')
    match = re.search(r'lastal', result)
    assert match is not None

    # the histogram tool from py-metagenomics
    result = get_version('histograms.py')
    match = re.search(r'py-metagenomics', result)
    assert match is not None

def test_apply_defaults():
    """ updating dicts without overwriting existing """

    # one level deep
    config = {"key1": "val1"}
    defaults = {"key1": "new_val1", "key2": "new_val2"}
    apply_defaults(config, defaults)
    assert config['key1'] == 'val1'
    assert config['key2'] == 'new_val2'

    # complicated nesting
    config = {'list1': [1, 2, 3],
              'dict1': {'key1': 'val1'},
              'dict2': {'key2': 'val2'},
              'dict3': {'key3': 'val3'},
             }
    defaults = {'list1': [4, 5],
                'dict1': {'key1': 'def_val1'},
                'dict2': {'key2b': 'def_val2'},
                'dict4': {'key4': 'def_val4'},
               }
    apply_defaults(config, defaults)
    assert config == {'list1': [1, 2, 3],
                      'dict1': {'key1': 'val1'},
                      'dict2': {'key2': 'val2',
                                'key2b': 'def_val2'},
                      'dict3': {'key3': 'val3'},
                      'dict4': {'key4': 'def_val4'},
                     }
