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
    for N in [1,2]:
        barcode = barcodes[N-1] if len(barcodes)>=N else barcodes[0]
        barcode_data["barcode{N}".format(N=N)] = barcode
        barcode_data["edocrab{N}".format(N=N)] = Seq.reverse_complement(barcode)
    return barcode_data

