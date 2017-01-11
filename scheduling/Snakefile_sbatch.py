#!/usr/bin/env python
"""
Submit this clustering script for sbatch to snakemake with:

    snakemake -j 99 --immediate-submit --cluster 'Snakefile-sbatch.py -t {threads} {dependencies}' --cluster-config cluster.yaml

The cluster config file should have values for sbatch options keyed byb snakemake rule name. A "default" rule cofig can be specified. EG:

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
                 ):
        self.scriptname = snakebashfile
        job_properties = read_job_properties(snakebashfile)
        print(job_properties, file=sys.stderr)
        self.rule = job_properties['rule']
        self.ifiles = job_properties['input']
        self.ofiles = job_properties['output']
        self.log_file = job_properties['params'].get('slurm_log',None)
        self.config = job_properties['cluster']
        self.threads = threads
        if dependencies == None or len(dependencies) < 1:
            self.dependencies = None
        else:
            # expects snakemake like list of numbers
            self.dependencies = dependencies
            assert len(self.dependencies) >= 1

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
                 ):
        SnakeJob.__init__(self, snakebashfile, 
                            threads=threads, 
                            dependencies=dependencies)

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
                    self.log_file=outfile
            else: 
                # a reasonable default for other rules
                self.log_file='logs/snakemake-{0}-slurm.out'.format(self.rule)
        else:
            # add slurm suffix so we don't collide with snakemake's use of log
            self.log_file+="-slurm.out"
        print("logging to {}".format(self.log_file), file=sys.stderr)
        # directory must exist, though:
        make_dir(os.path.dirname(os.path.abspath(self.log_file)))

        # options from config file
        configured_options=""
        # print(str(self.config), file=sys.stderr)
        for key, value in self.config.items():
            if key=='freeform':
                configured_options += " " + value
            else:
                flag = key
                configured_options += " " + flag + " " + str(value)
                # print("Adding {}:{}".format(key,value),file=sys.stderr)

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dependencies", nargs="*", \
                        help="{{dependencies}} string given by snakemake\n")
    parser.add_argument("-t","--threads", type=int, default=1, \
                        help="{{threads}} string given by snakemake\n")
    parser.add_argument("snakescript", help="Snakemake generated shell \
                         script with commands to execute snakemake rule\n")
    args = parser.parse_args()

    #print("Passed bidniz:", args.snakescript, args.dependencies, file=sys.stderr)
    #print("Passed args:", args, file=sys.stderr)
    sj = SnakeJobSbatch(args.snakescript, \
                        threads=args.threads,\
                        dependencies=args.dependencies, 
                        )
    try:
        sj.schedule()
    except UndefinedJobRule as err:
        print(err.msg, file=sys.stderr)
        sys.exit(2)
