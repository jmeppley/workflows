import os
import subprocess

def get_bbduk_resource_path(config):
    """
    find the bbduk.sh script and look for the rsources folder
    """

    # if there is a bbmap path set, look there
    if 'bbmap' in config.get('system_paths',{}):
        resource_path = config['system_paths']['bbmap'] + "/resources"
        if os.path.exists(resource_path) and os.path.isdir(resource_path):
            return resource_path

    # infer path from binary location (using which and subprocess)
    try:
        cmd='find `which bbduk.sh | sed s/bbduk\.sh/../` ' + \
             '-name "adapters.fa" | sed s#/adapters.fa##'
        resource_path = \
            subprocess.check_output(cmd, shell=True).decode().strip()
    except:
        print("Can't find bbmap resources, please set system_paths->bbmap in the snakemake config file")
        raise
    else:
        return resource_path

