#!/usr/bin/env python
"""
Submit this clustering script for sbatch to snakemake with:

    snakemake -j 99 --debug --immediate-submit --cluster 'Snakefile-sbatch.py  -t {threads} {dependencies}'
"""
import argparse
import sys
import subprocess
import os
import math
import errno
import json
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
                 config=None):
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
        self.config = config

class UndefinedJobRule(Exception):
    """Exception in case an sbatch job has no defined resource usage in the
    code."""
    def __init__(self, msg):
        self.msg = msg


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
                self.log_file=self.ofiles[0] + '-slurm.out'
            else: 
                # a reasonable default for other rules
                self.log_file='snakemake-{0}-slurm.out'.format(self.rule)
        # directory must exist, though:
        make_dir(os.path.dirname(os.path.abspath(self.log_file)))

        schedule_rule = "schedule_{0}".format(self.rule)
        if schedule_rule in self.config:
            rule_conf = self.config[schedule_rule]
            # If rule_conf is referring to another scheduling rule, use those
            # resources instead
            try:
                if rule_conf.startswith("schedule_"):
                    rule_conf = self.config[rule_conf]
            except KeyError:
                raise UndefinedJobRule('No schedule config found for {0}'.format(rule_conf))
            except AttributeError:
                pass

            

            attributes = {
                    'dep_str': self.dep_str,
                    'job_name': 'snakemake_{0}'.format(self.rule),
                    'sbatch_job_path': self.config['sbatch_general']['wrapper_script'],
                    'script_name': self.scriptname,
                    'days': rule_conf['days'],
                    'hours': rule_conf['hours'],
                    'minutes': rule_conf['minutes'],
                    'partition': rule_conf['partition'],
                    'cores': rule_conf.get('cores',self.threads),
                    'account': self.config['sbatch_general']['account'],
                    'log_file': self.log_file,
                    'extra_parameters': rule_conf.get('extra_parameters', "")
            }
            sbatch_cmd = """sbatch --output={log_file} {dep_str} -A {account} -p {partition} -n {cores} -t {days}-{hours}:{minutes}:00 \
                    '{script_name}'""".format(**attributes)
        else:
            raise UndefinedJobRule('No schedule config found for schedule_{0}'.format(self.rule))
            return 2

        print(sbatch_cmd, file=sys.stderr)
        popenrv = subprocess.Popen(sbatch_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).communicate()

        # Snakemake expects only id of submitted job on stdout for scheduling
        # with {dependencies}
        try:
            print("%i" % int(popenrv[0].split()[-1]))
        except ValueError:
            print("Not a submitted job: %s" % popenrv[0])
            sys.exit(2)

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
    sj = SnakeJobSbatch(args.snakescript, \
                        threads=args.threads,\
                        dependencies=args.dependencies, \
                        config=json.load(open(args.configfile)))
    try:
        sj.schedule()
    except UndefinedJobRule as err:
        print(err.msg, file=sys.stderr)
        sys.exit(2)
