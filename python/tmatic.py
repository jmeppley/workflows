"""
methods and constants to help with runnning trimmomatic
"""
from Bio import Seq

primer_templates = {
    'truseq':""">PrefixTruA/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixTruA/2
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC{edocrab2}ATCTCGTATGCCGTCTTCTGCTTG
>PrefixTruB/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixTruB/2
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC{edocrab2}CGATCTCGTATGCCGTCTTCTGCTTG""",
    'scripseq':""">PrefixScript/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixScript/2
CAAGCAGAAGACGGCATACGAGAT{edocrab2}GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT""",
    'nextera':""">Prefix25/1
AATGATACGGCGACCACCGAGATCTACAC{barcode1}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Prefix25/2
CAAGCAGAAGACGGCATACGAGAT{edocrab2}GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG""",
}

def get_primer_template(chemistry):
    return primer_templates[chemistry]

def process_barcode_string(barcode_string):
    barcodes = barcode_string.split(":")
    barcode_data = {}
    for N in [1, 2]:
        barcode = barcodes[N-1] if len(barcodes)>=N else barcodes[0]
        barcode_data["barcode{N}".format(N=N)] = barcode
        barcode_data["edocrab{N}".format(N=N)] = Seq.reverse_complement(barcode)
    return barcode_data

def get_chemistry_barcodes(sample, config):
    """
    Attempt to work out cleaning params for given sample.

    Try the following:
        * look in config[sample_data][{sample}] for
          * sample_sheet
          * chemistry and barcodes
        * use Ns for barcode and use config.get(chemistry, 'scripseq')
    """
    sample_data = config['sample_data'][sample]
    if 'barcode' in sample_data and 'chemistry' in sample_data:
        return sample_data['chemistry'], [sample_data['barcode'], ]
    if 'barcodes' in sample_data and 'chemistry' in sample_data:
        return sample_data['chemistry'], sample_data['barcodes']
    if 'sample_sheet' in sample_data:
        return parse_sample_sheet(sample, sample_data['sample_sheet'])
    else:
        if config.get('chemistry', 'scripseq') == 'nextera':
            return config['chemistry'], ['NNNNNN', 'NNNNNN']
        return config.get('chemistry', 'scripseq'), ['NNNNNN', ]

def parse_sample_sheet(sample, sample_sheet):
    """
    return chemistry and barcode for this sample
    """
    #TODO
    raise NotImplementedError("You'll have to set the barcode and chemsitrya"
                              "manually for now!")

