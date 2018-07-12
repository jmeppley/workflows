######
# collection of methods for processing riboframe output
#####

import os, sys, re
from Bio import SeqIO

region_map={'V1': (69, 99),
         'V2': (137, 242),
         'V3': (433, 497),
         'V4': (576, 682),
         'V5': (822, 879),
         'V6': (986, 1043),
         'V7': (1117, 1173),
         'V8': (1243, 1294),
         'V9': (1435, 1465)}

def get_region(region_string, region_map=region_map):
    if region_string in region_map:
        return region_map[region_string]
    else:
        try:
            start, end = re.search(r'(\d+)[^0-9]+(\d+)', region_string).groups()
            #Jess B changed the code here:
            return (int(region_map['V'+start][0]), int(region_map['V'+end][1]))
        except:
            raise Exception("Cannot guess start/end from region string: %s" % (region_string))

def get_ribotrap_16S_region_sequence(record, 
                                     padding=0,
                                     region_map=region_map):
    """
    Parse the record ID to get the 16S coordinates, and return copy of record clipped to the matching region
    
    The record id should end in something like ...1234:123.189 (where 123 and 189 are the start and end positions)
    """
    start,end = sorted([int(p) for p in record.id.split(":")[-1].split('.')])
    for region, (rstart,rend) in region_map.items():
        if rstart>=start and rend<=end:
            # Region MATCH!
            break
    else:
        return None, None
    
    # clip record
    clip_start = max(rstart-start-padding,0)
    clip_end = min(rend-start+padding,len(record))
    return region,record[clip_start:clip_end]

def pull_out_region_piped(
                     region,
                     handle=None,
                     region_map=region_map,
                     renaming_root=None,
					 padding=0,
                    ):

    # initialize counters
    region_count=0
    
    # initialize output files (one for each region)
    if handle is None:
        handle=sys.stdout
        
    # loop over input records
    for record in SeqIO.parse(sys.stdin, 'fasta'):
        r, subrecord = get_ribotrap_16S_region_sequence(
                record, 
                region_map={region: get_region(region, region_map),                           },
                padding=padding)
        if r is not None:
            region_count+=1
            if renaming_root is not None:
                subrecord.id = "%s_%d" % (renaming_root, region_count)
            SeqIO.write([subrecord,],handle,'fasta')

    return

def pull_out_regions(ribotrap_fasta, 
                     regions=['V6','V3'],
                     outfile_template=None,
                     region_map=region_map,
                     renaming_root=None,
                    ):

    # initialize counters
    region_counts = {r:0 for r in regions}
    region_counts[None]=0
    
    # initialize output files (one for each region)
    if outfile_template is None:
        file_name, file_ext = os.path.splitext(ribotrap_fasta)
        outfile_template =  file_name + ".%s" + file_ext
    for region in regions:
        outfile_name = outfile_template % (region)
        handles[region]=open(outfile_name,'w')
        
    # loop over input records
    for record in SeqIO.parse(ribotrap_fasta, 'fasta'):
        region, subrecord = get_ribotrap_16S_region_sequence(
                record, 
                region_map={r:get_region(r, region_map) for r in regions},
                )
        region_counts[region]+=1
        if region is not None:
            if renaming_root is not None:
                subrecord.id = "%s_%d" % (renaming_root, region_counts[region])
            SeqIO.write([subrecord,],handles[region],'fasta')

    # close output handles
    for handle in handles.values():
        handle.close()
        
    # return counts
    return region_counts


