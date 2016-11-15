#!/usr/bin/env python
"""
Submit this clustering script for sbatch to snakemake with:

    snakemake -j 99 --debug --immediate-submit --cluster 'Snakefile-sbatch.py [-c CONFIG_FILE}] -t {threads} {dependencies}'

The sbatch config file should have values for sbatch options keyed byb snakemake rule name. A "default" rule cofig can be specified. EG:

default:
    --time: 0-5:0
    --cpus-per-task: 1
    --mem-per-cpu: 1G
    --partition: scope-shared.q
heavy-rule:
    --time: 14-0
    --partition: lm.q
    --cpus-per-task: 40
    --mem-per-cpu: 0
ram-only-rule:
    --mem-per-cpu: 100G

Any rule (or the defualt) can have a "freeform" key with options to pass

"""
import argparse
import sys
import subprocess
import os
import math
import errno
import json
import yaml
import re
from snakemake.utils import read_job_properties

def make_dir(directory):
    """Make directory unless existing. Ignore error in the latter case."""
    try:
        os.makedirs(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


class SnakeJob:
    """Snakemake can generate bash scripts that can be sumbitted by a
    scheduler.  This class reads the bash script and stores the number of the
    rule, name of bash file and the supplied input files."""
    def __init__(self, snakebashfile, \
                 threads=1, \
                 dependencies=None, \
                 config={}):
        self.scriptname = snakebashfile
        job_properties = read_job_properties(snakebashfile)
        self.rule = job_properties['rule']
        self.ifiles = job_properties['input']
        self.ofiles = job_properties['output']
        self.log_file = job_properties['params'].get('slurm_log',None)
        self.threads = threads
        if dependencies == None or len(dependencies) < 1:
            self.dependencies = None
        else:
            # expects snakemake like list of numbers
            self.dependencies = dependencies
            assert len(self.dependencies) >= 1

        # get config for this rule, falling back to any defaults
        self.config = config.get('default',{})
        self.config.update(config.get(self.rule,{}))

class UndefinedSbatchFlag(Exception):
    """Exception in case a rule is configred to use a flag that doesn't exist"""
    def __init__(self, msg):
        self.msg = msg

class UndefinedJobRule(Exception):
    """Exception in case an sbatch job has no defined resource usage in the
    code."""
    def __init__(self, msg):
        self.msg = msg

sbatch_flags=set((
        '--account',
        '--array',
        '--bb',
        '--begin',
        '--clusters',
        '--comment',
        '--constraint',
        '--core-spec',
        '--cores-per-socket',
        '--cpu-freq',
        '--deadline',
        '--dependency',
        '--distribution',
        '--error',
        '--exclude',
        '--export-file',
        '--extra-node-info',
        '--gid',
        '--gres',
        '--gres-flags',
        '--hint',
        '--input',
        '--jobid',
        '--job-name',
        '--licenses',
        '--mail-type',
        '--mail-user',
        '--mcs-label',
        '--mem',
        '--mem_bind',
        '--mem-per-cpu',
        '--mincpus',
        '--nodefile',
        '--nodelist',
        '--nodes',
        '--ntasks',
        '--ntasks-per-core',
        '--ntasks-per-node',
        '--ntasks-per-socket',
        '--output',
        '--partition',
        '--power',
        '--priority',
        '--profile',
        '--qos',
        '--reservation',
        '--signal',
        '--sockets-per-node',
        '--switches',
        '--thread-spec',
        '--threads-per-core',
        '--time',
        '--time-min',
        '--tmp',
        '--uid',
        '--wckey',
        '--workdir',
))

class SnakeJobSbatch(SnakeJob):
    def __init__(self, snakebashfile, \
                 threads=1, \
                 dependencies=None, \
                 config=None):
        SnakeJob.__init__(self, snakebashfile, threads, dependencies, config)

        # Set up dependency string for SLURM
        if self.dependencies == None:
            self.dep_str = ''
        else:
            self.dep_str = '-d ' + ','.join(["afterok:%s" % d for d in self.dependencies])

    def schedule(self):
        """Schedules a snakemake job with sbatch and determines resource usage
        based on input files."""
        # set up SLRUM log file
        if self.log_file is None:
            if len(self.ofiles)>0:
                # otherwise use the first output file
                outfile=self.ofiles[0] +'-slurm.out'
                outdir,outname=os.path.split(outfile)
                if len(outdir.strip())==0:
                    self.log_file=os.path.join('logs',outname)
                else:
                    self.log_file=os.path.join(outdir,'logs',outname)
            else: 
                # a reasonable default for other rules
                self.log_file='logs/snakemake-{0}-slurm.out'.format(self.rule)
        else:
            # add slurm suffix so we don't collide with snakemake's use of log
            self.log_file+="-slurm.out"
        print(self.log_file)
        # directory must exist, though:
        make_dir(os.path.dirname(os.path.abspath(self.log_file)))

        # options from config file
        configured_options=""
        for key, value in self.config.items():
            if key=='freeform':
                configured_options += " " + value
            else:
                flag = key
                # validate flags?
                #if flag not in sbatch_flags:
                #    raise UndefinedSbatchFlag("{} is not a recognized flag for sbatch")
                #    return 2
                configured_options += " " + flag + " " + value

        attributes = {
                'dep_str': self.dep_str,
                'job_name': 'snakemake_{0}'.format(self.rule),
                'script_name': self.scriptname,
                'log_file': self.log_file,
                'configured_options': configured_options,
        }
        sbatch_cmd = """sbatch --open-mode=append --output={log_file} \
                -J {job_name} {dep_str} {configured_options} \
                '{script_name}'""".format(**attributes)

        print(sbatch_cmd, file=sys.stderr)
        popenrv = subprocess.Popen(sbatch_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).communicate()

        # Snakemake expects only id of submitted job on stdout for scheduling
        # with {dependencies}
        try:
            slurm_job_id=int(popenrv[0].split()[-1])
            print("%i" % slurm_job_id)
        except ValueError:
            print("Not a submitted job: %s" % popenrv[0])
            sys.exit(2)

        # write JOB ID to file for easier tracking
        with open(self.log_file + ".jid",'w') as JIDOUT:
            JIDOUT.write(str(slurm_job_id))
            JIDOUT.write('\n')

def load_config(configfile):
    with open(configfile) as CONFIG:
        if os.path.splitext(configfile)[-1]=='.yaml':
            return yaml.load(CONFIG)
        else:
            return json.load(CONFIG)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dependencies", nargs="*", \
                        help="{{dependencies}} string given by snakemake\n")
    parser.add_argument("-c","--configfile", default="config_sbatch.json", \
                        help="{{threads}} string given by snakemake\n")
    parser.add_argument("-t","--threads", type=int, default=1, \
                        help="{{threads}} string given by snakemake\n")
    parser.add_argument("snakescript", help="Snakemake generated shell \
                         script with commands to execute snakemake rule\n")
    args = parser.parse_args()

    #print("Passed bidniz:", args.snakescript, args.dependencies, file=sys.stderr)
    #print("Passed args:", args, file=sys.stderr)
    config=load_config(args.configfile)
    sj = SnakeJobSbatch(args.snakescript, \
                        threads=args.threads,\
                        dependencies=args.dependencies, \
                        config=config)
    try:
        sj.schedule()
    except UndefinedJobRule as err:
        print(err.msg, file=sys.stderr)
        sys.exit(2)
