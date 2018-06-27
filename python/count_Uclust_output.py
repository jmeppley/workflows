# 
#
# Copyright (C) OPEN ACCESS
# Author: Jessica Bryant
# URL: <https://github.com/JessAwBryant/>

"""
Tabling Uclust output

"""

import re as re
import sys as sys
from collections import defaultdict
import pandas as pd

D3_level = defaultdict(lambda: 0, {})
D5_level = defaultdict(lambda: 0, {})
D6_level = defaultdict(lambda: 0, {})

input_file = sys.argv[1] 

if re.search('18S', input_file):
    domain = 'Eukaryote'
    
if re.search('16S', input_file):
    domain = 'Bacteria'
        
output_file_name1 = input_file.strip(".txt") + ".D3.counts.txt"
output_file_name2 = input_file.strip(".txt") + ".D5.counts.txt"
output_file_name3 = input_file.strip(".txt") + ".D6.counts.txt"

f = open(input_file, 'r')
#g = open(output_file_name, 'r')

for input_line in f:
    line = input_line.split('\t')[1].split(';')
    line = [re.sub(r'D_[0-9].*__', '', x) for x in line]
    line = [x for x in line if not re.search(x, ' +')]
    
    line = line + ['Unassigned'] * (7-len(line))
    
    if line[0] == 'Eukaryota' and domain == 'Bacteria':
        continue
    
    if line[0] == 'Bacteria' and domain == 'Eukaryote':
        continue

    D3 = line[3]
    D5 = line[5]
    D6 = line[6]

    if D5 == '':
        print(input_line)
    
    #make some tweeks to D5
    if re.search('uncultured|Ambiguous|Unknown|unknown|Unassigned', D5):
        D5 = line[4]

    if re.search('uncultured|Ambiguous|Unknown|unknown|Unassigned', D5):
        D5 = line[3]

    if re.search('uncultured|Ambiguous|Unknown|unknown|Unassigned', D5):
        D5 = line[2]

    if re.search('uncultured|Ambiguous|Unknown|unknown|Unassigned', D5):
        D5 = line[1]

    # if the ID is to SAR202
    if line[3] == 'SAR11 clade':
        D3 = 'SAR11 clade'
        D5 = 'SAR11 clade'
            
    # if the ID is to SAR202
    if line[2] == 'SAR202 clade':
        D3 = 'SAR202 clade'
        D5 = 'SAR202 clade'
        D6 = 'SAR202 clade'

    if line[2] == 'Chloroplast':
        D3 = 'Chloroplast'
        D5 = 'Chloroplast'
        D6 = 'Chloroplast'

    if D5 == 'uncultured thaumarchaeote':
       D5 = 'Marine Group I'
    
    if line[1] == 'Marinimicrobia (SAR406 clade)':
        D3 = 'Marinimicrobia'
        D5 = 'Marinimicrobia'
        D6 = 'Marinimicrobia'

    if D5 == 'SAR324 clade(Marine group B)':
        D5 = 'SAR324 clade'
        D6 = 'SAR324 clade'

    if D5 == 'Prochlorococcus':
        D5 = 'FamilyI'

    D1 = line[1]
    if D1 == 'Thaumarchaeota':
        D5 = 'Thaumarchaeota'
    
    D3_level[D3] = D3_level[D3]+1
    D5_level[D5] = D5_level[D5]+1
    D6_level[D6] = D6_level[D6]+1
    
pd.DataFrame(D3_level, index=[0]).transpose().to_csv(output_file_name1, sep='\t')
pd.DataFrame(D5_level, index=[0]).transpose().to_csv(output_file_name2, sep='\t')
pd.DataFrame(D6_level, index=[0]).transpose().to_csv(output_file_name3, sep='\t')
