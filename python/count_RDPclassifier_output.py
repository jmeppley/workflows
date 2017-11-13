# 
#
# Copyright (C) OPEN ACCESS
# Author: Jessica Bryant
# URL: <https://github.com/JessAwBryant/>

"""
Tabling RDP output
"""

import re as re
import sys as sys
from collections import defaultdict
import pandas as pd

D3_level = defaultdict(lambda: 0, {})
D6_level = defaultdict(lambda: 0, {})

input_file = sys.argv[1] 
output_file_name = input_file.strip(".txt") + ".counts.txt"
f = open(input_file, 'r')
g = open(output_file_name, 'r')

for input_line in f:
    line=input_line.split('\t')[1].split(';')
    line = line + ['unknown'] * (7-len(line))
    
    D3 = line[3].strip("D_").strip("3__")
    D6 = line[6].strip("D_").strip("6__")
    
    D3_level[D3] = D3_level[D3]+1
    D6_level[D6] = D6_level[D6]+1
    
pd.DataFrame(D3_level, index=[0]).transpose().to_csv()
pd.DataFrame(D6_level, index=[0]).transpose().to_csv(output_file_name, sep='\t')
