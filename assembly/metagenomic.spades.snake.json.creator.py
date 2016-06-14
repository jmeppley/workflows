#!/usr/bin/env python3
"""
Take a json file and replace the inputs section with a custom set

This works two ways:

    (1) give a list of runs and linst of sample nums and create a single config file using all samples
    (2) give a list of runs and create a config file for each unique sample name in those runs
"""

import argparse, glob, os, sys, re, json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--run", action="append", help="Run name(s) to pull reads from")
    parser.add_argument("-s", "--sample_num", action="append", help="Sample num to pull from run. Specify once for each run name given. If you desire mutiple samples from a given run, list the run multiple times")
    parser.add_argument("-R", "--run_root", default=".", help="Folder in which to find run names")
    parser.add_argument("-i", "--input_json", help="JSON file to modify")
    parser.add_argument("-o", "--output_json", help="JSON file to write")
    args = parser.parse_args()

    # switch to absolute path
    args.run_root = os.path.abspath(args.run_root)

    # LOAD template
    if args.input_json is not None:
        with open(args.input_json) as INF:
            config = json.load(INF)
    else:
        config = json.load(sys.stdin)

    # CLEAR dummy input data
    config["inputs"]={}
    
    # find read files and update data structure
    if args.sample_num is not None:
        # create a single JSON file
        for run_name, sample_num in zip(args.run, args.sample_num):
            sample_data = {"sample_num":sample_num}
            run_dir = os.path.join(args.run_root, run_name)
            for root, dirs, files in os.walk(run_dir):
                for file_name in files:
                    if file_name=='SampleSheet.csv':
                        sample_data['sample_sheet']=os.path.join(root,file_name)
                        continue
                    m = re.search(r'^([^_]+)_S(\d+)_.*R([12])_',file_name)
                    if m:
                        sample_name=m.group(1)
                        sample_num=int(m.group(2))
                        dir = "R1" if m.group(3)=="1" else "R2"
                        sample_data[dir]=os.path.join(root,file_name)
                        continue
                if sample_name is None:
                    raise Exception("No fastq found for sample %d in %s" % (sample_num, run_name))
                config['inputs'][sample_name]=sample_data

        # OPEN output handle
        if args.output_json is not None:
            OUTF = open(args.output_json,'w')
        else:
            OUTF = sys.stdout

        # WRITE JSON
        json.dump(config, OUTF, indent=1, sort_keys=True)

        # CLOSE output if needed
        if args.output_json is not None:
            OUTF.close()

    else:
        # find all unique sample names
        sample_data_dict={}
        sample_sheet_dict={}
        for run_name in args.run:
            run_dir = os.path.join(args.run_root, run_name)
            for root, dirs, files in os.walk(run_dir):
                for file_name in files:
                    if file_name=='SampleSheet.csv':
                        sample_sheet_dict[run_name]=os.path.join(root,file_name)
                        continue
                    m = re.search(r'^([^_]+)_S(\d+)_.*R([12])_',file_name)
                    if m:
                        sample_name=m.group(1)
                        sample_num=int(m.group(2))
                        dir = "R1" if m.group(3)=="1" else "R2"
                        sample_data_dict.setdefault(sample_name,{}).setdefault(run_name,{}).setdefault(sample_num,{})[dir]=os.path.join(root,file_name)
                        continue

        if args.output_json is not None:
            out_root, out_name = os.path.split(args.output_json)
        else:
            if args.input_json is not None:
                out_root, out_name = os.path.split(args.input_json)
            else:
                raise Exception("CAnnot write multiple JSON files to STDOUT")

        for (sample_name, run_data) in sample_data_dict.items():
            sample_count=0
            config['inputs']={}
            for run_name, sample_num_data in run_data.items():
                if run_name not in sample_sheet_dict:
                    raise Exception("No sample sheet found for run %s" % run_name)
                sample_sheet = sample_sheet_dict[run_name]
                for sample_num, file_data in sample_num_data.items():
                    sample_count+=1
                    if sample_count>1:
                        sample_key="%s_%d" % (sample_name,sample_count)
                        if sample_count==2:
                            #rename first one
                            config['inputs']["%s_1"%(sample_name)] = \
                                    config['inputs'][sample_name]
                            del config['inputs'][sample_name]
                    else:
                        sample_key=sample_name
                    file_data['sample_num']=sample_num
                    file_data['sample_sheet']=sample_sheet
                    config['inputs'][sample_key]=file_data

            # OPEN output handle
            if not os.path.exists(os.path.join(out_root, sample_name)):
                os.makedirs(os.path.join(out_root, sample_name))
            out_file=os.path.join(out_root, sample_name,out_name)

            # WRITE JSON
            with open(out_file,'w') as OUTF:
                json.dump(config, OUTF, indent=1, sort_keys=True)

if __name__ == "__main__":
    main()
