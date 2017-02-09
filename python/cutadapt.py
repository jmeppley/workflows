from Bio import Seq

# The scripseq adapters determined by trial and error
#  The others are best guesses (revcomp of end of /1 and start of /2)
adapters = {
        'truseq': ['GCTCTTCCGATC', 'AGATCGGAAGAG'].
        'scripseq': ['GATCGGAAGAGC', 'GTCTTCTGCTTG'],
        'nextera': ['CTGTCTCTTATA', 'GTCTTCTGCTTG'],
        }

def get_adapters(chemistry):
    return adapters[chemistry]
