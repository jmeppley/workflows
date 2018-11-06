" Tests for python/common.py "
import os
import python.common

def test_in_working_dir():
    " should return for files in repo "
    # files in repo
    test_path = "test/bats"
    assert python.common.is_in_working_dir(test_path)

    test_path = os.path.abspath("python/common.py")
    assert python.common.is_in_working_dir(test_path)

    repo_dir = os.path.basename(os.getcwd())
    test_path = "../{}/common/gunzip.snake" \
                                .format(repo_dir)
    assert python.common.is_in_working_dir(test_path)

    test_path = os.path.abspath(test_path)
    assert python.common.is_in_working_dir(test_path)

    #files not in repo
    for test_path in ["/bin/ls", "/home", "/Users"]:
        if os.path.exists(test_path):
            assert not python.common.is_in_working_dir(test_path)

def test_get_version():
    " try some things expected to be in path "
    assert python.common.get_version('python')\
                .lower()\
                .startswith('python')

def test_parse_stats():
    test_file = "test/bats/outputs/nose/fortyk.R1.fastq.stats"
    assert (python.common.parse_stats(test_file) ==
            {'reads': 10000, 'bases': 1510000})

def test_apply_defaults():
    " make sure this works as intended "
    config = {
        "param_1": True,
        "param_2": False,
        "subset": {
            "param_1": True,
            "param_2": False,
        }
    }

    defaults = {"param_2": True, "param_3": True}
    python.common.apply_defaults(config, defaults)
    assert (config ==
            {
                "param_1": True,
                "param_2": False,
                "param_3": True,
                "subset": {
                    "param_1": True,
                    "param_2": False,
                }
            })

    defaults = {"subset":{"param_2": True, "param_3": True}}
    python.common.apply_defaults(config, defaults)
    assert (config ==
            {
                "param_1": True,
                "param_2": False,
                "param_3": True,
                "subset": {
                    "param_1": True,
                    "param_2": False,
                    "param_3": True,
                }
            })

def test_get_file_name():
    " simple utility to mak sure we get a string "
    file_name = "/some/file/path"
    assert python.common.get_file_name(file_name) == file_name
    assert python.common.get_file_name([file_name,]) == file_name
    assert python.common.get_file_name((file_name,)) == file_name

def test_get_set_from_file():
    " turns file into set() of trimmed lines "
    file_name = "test/conda/assembly.yml"
    data = python.common.get_set_from_file(file_name)
    assert isinstance(data, set)
    assert "channels:" in data
    assert "- bioconda" in data
    assert "foobar" not in data
