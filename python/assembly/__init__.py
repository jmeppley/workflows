"""
Functions used in the assembly workflows
"""
from python.assembly.stats import \
        get_contig_stats, \
        generate_histogram, \
        get_contig_length_summary_stats

from python.assembly.coverage import \
        get_annot_coverage_stats, \
        get_coverage_stats

from python.assembly.annotation import \
        filter_annotations, \
        drop_all_rna_overlaps, \
        filter_and_extract_rRNA


"""
## I think this is new code that I forgot to finish...
import pandas, re
def orphan_code(annot_file):
    "" What was this supposed to do?
    It turnd a gff file into a table of gene positions...

    I think it' been replaced by annotations.get_gene_name
    ""

    genes = pandas.read_table(annot_file,
                              usecols=[0, 1, 2, 3, 4],
                              names=['contig', 'program', 'model',
                                     'start', 'end'])

    gene_col = []
    current_contig = None
    for contig, program, model in genes[['contig',
                                         'program',
                                         'model']].values:
        if contig != current_contig:
            gene_count, rrna_count, trna_count = 0, 0, 0
            current_contig = contig

        # is it a CDS or RNA
        if program.lower().startswith('prodigal') and model == 'CDS':
            gene_count += 1
            suffix = "_{}".format(gene_count)
        elif program == 'cmsearch' and re.search('rRNA', model) is not None:
            rrna_count += 1
            suffix = "_rRNA_{}".format(rrna_count)
        elif program == 'cmsearch' and re.search('tRNA', model) is not None:
            trna_count += 1
            suffix = "_tRNA_{}".format(trna_count)
        else:
            raise Exception("Can't classify GFF line:\n" +
                            repr((contig, program, model)))

        gene_col.append(contig + suffix)

    genes['gene'] = gene_col
    return genes[['gene', 'contig', 'start', 'end']]
"""
