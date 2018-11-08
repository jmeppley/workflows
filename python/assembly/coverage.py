"""
Functions for calculating coverages in assembly workflows:
    get_annot_coverage_stats: make table of covs for each annot
    get_coverage_stats: make table of covs for each contig

Mostly get data from .depths files
"""
from collections import defaultdict
import numpy
import pandas
from python.assembly.annotation import get_mixed_annotations
from Bio import SeqIO
from snakemake import logger

def get_coverage_stats(contig_depth_file,
                       contig_fasta,
                       contig_read_counts_file,
                       contig_stats_out):
    """
    Generate table of read counts and mapped coverage from samtools output
    """
    print("getting coverage stats")
    # add other files if requested
    # read counts
    logger.info("Parsing read count file: {}"
                .format(contig_read_counts_file))
    read_count_table = pandas.read_table(contig_read_counts_file,
                                         delim_whitespace=True,
                                         names=['ReadCount',
                                                'Contig'])\
                             .set_index('Contig')

    # convert base by base depth data into coverage
    logger.info("Parsing read depth file: {}"
                .format(contig_depth_file))
    mapping_depth_table = get_samtool_depth_table(contig_depth_file,
                                                  contig_fasta,
                                                 )
    contig_stats = mapping_depth_table.join(read_count_table,
                                            how='left').fillna(0)

    for col in ['Length', 'ReadCount',
                'MaxCov', 'MinCov', 'CumuLength']:
        if col in contig_stats.columns:
            contig_stats[col] = contig_stats[col].astype(int)

    logger.info("Writing coverage table to: %s"
                 % (contig_stats_out))
    contig_stats.to_csv(contig_stats_out,
                        sep='\t',
                        float_format="%0.2f")


def get_annot_coverage_stats(contig_depth_file,
                             annot_gff_file,
                             annot_cov_out):
    """
    Generate table of mapped coverage for annotations
    """
    contig_annots = get_mixed_annotations(annot_gff_file)

    with open(annot_cov_out, 'wt') as out_handle:
        out_handle.write("gene\tcoverage\n")
        with open(contig_depth_file) as depth_handle:
            contig_depths = contig_depths_generator(depth_handle)
            for contig, depths in contig_depths:
                try:
                    contig_genes = contig_annots[contig]
                except KeyError:
                    # no genes for contig
                    continue
                contig_extent = max(max(g[1:]) \
                                    for g in contig_genes)
                coverage = \
                        numpy.array(list(_insert_zeros(depths,
                                                       contig_extent)))
                for gene, start, end in contig_genes:
                    start, end = sorted([start, end])
                    gene_cov = coverage[start - 1:end].mean()
                    out_handle.write("{}\t{:0.2f}\n"
                                     .format(gene, gene_cov))


def get_samtool_depth_table(depth_file, fasta_file=None):
    """
    Calculate coverage stats for each contig in an assembly

    Params:
     depth_file: output file from the command:
                    `samtools depth reads.v.contigs.bam`
                 this is a 3 column file with one line per base.
                 columns are:
                     'contig_id base_index base_depth'
    fasta_file: fasta of contigs, used to get length (optional)

    Returns:
     pandas.DataFrame with one row per contig and the following columns:
            Contig, MeanCov, MedCov, MaxCov, MinCov, StdCov
            """
    contig_lengths = {r.id:len(r) for r in
                      SeqIO.parse(fasta_file, 'fasta')} \
                     if fasta_file is not None \
                     else None

    with open(depth_file) as depth_handle:
        return get_contig_coverage_table(depth_handle, contig_lengths)


CONTIG_COVERAGE_COLUMNS = ['Contig',
                           'MedCov',
                           'Q2Q3Cov',
                           'MeanCov',
                           'StdCov',
                           'MinCov',
                           'MaxCov',
                          ]

def get_contig_coverage_table(contig_depths_handle,
                              contig_lengths=None):
    """
    Uses this parser to generate a DataFrame of coverage stats by contig
    Pass in dict of contig legths if 0 depths are excluded from depth file
    """
    if contig_lengths is None:
        contig_lengths = {}

    table = \
        pandas.DataFrame(_contig_cov_row_generator(contig_depths_handle,
                                                   contig_lengths),
                         columns=CONTIG_COVERAGE_COLUMNS) \
              .set_index('Contig')
    return table


def _contig_cov_row_generator(contig_depths_handle, contig_lengths):
    depths_iterator = contig_depths_generator(contig_depths_handle)
    for contig, depths in depths_iterator:
        coverage = \
            numpy.array(list(_insert_zeros(depths,
                                           contig_lengths.pop(contig,
                                                              0))))
        yield (contig,
               numpy.median(coverage),
               q2q3_mean(coverage),
               coverage.mean(),
               coverage.std(),
               coverage.min(),
               coverage.max())

    for contig in contig_lengths:
        yield (contig, 0, 0, 0, 0, 0, 0)

def contig_depths_generator(depth_file_handle):
    """ yield a tuple contining contig and list of depth entries """
    last_contig = None
    for line in depth_file_handle:
        contig, depth_data = _parse_depth_line(line)

        # collect lines until we see a new contig
        if last_contig != contig:
            if last_contig is not None:
                yield last_contig, contig_depths
            last_contig = contig
            contig_depths = []

        contig_depths.append(depth_data)

    # yield final contig
    if last_contig is not None:
        yield last_contig, contig_depths

def _parse_depth_line(depth_line):
    """ return contig, (base, depth) nested tuples with base & depth as
    integers """
    data = depth_line.rstrip().split('\t')
    return data[0], tuple(int(d) for d in data[1:])

def _insert_zeros(contig_depths, contig_length=0):
    """
    The depths file only lists bases with non-zero depths

    return only depth values (not bases), but insert zeros where needed
    """

    # keep track of the last base with a >0 depth
    last_base = 0

    # loop over depth values
    for base, depth in contig_depths:
        last_base += 1
        while base > last_base:
            # insert 0's until we catch up
            yield 0
            last_base += 1
        yield depth

    # insert 0's until we reach the end
    while last_base < contig_length:
        last_base += 1
        yield 0

def q2q3_mean(coverage):
    """ Calculate the mean of the 2nd and 3rd quartiles """
    if len(coverage) < 4:
        return coverage.mean()
    coverage.sort()
    quartile_size = int(len(coverage) / 4)
    return coverage[quartile_size:-quartile_size].mean()
