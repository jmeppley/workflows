"""
Tests for cmsearch.py
"""
import re
from python.cmsearch import get_cm_model_version

def test_get_cm_model_version():
    """ get the version of the bundled models """
    version = get_cm_model_version('./models')
    assert re.search('Release', version) is not None

