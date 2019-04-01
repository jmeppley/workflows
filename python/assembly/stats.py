"""
Functions used to process stats in the assembly workflows"
    get_contig_stats: get GC and length from contigs
    generate_histogram:Generate ASCII histrogram of the selected data
    get_contig_length_summary_stats: N50, mean length, etc

"""
import numpy
import pandas
from Bio import SeqIO, SeqUtils
from edl.util import ascii_histogram
from snakemake import logger

def get_contig_stats(contigs_fasta,
                     contig_stats_out=None):
    """
    Extracts GC and length from contigs fasta
    """
    # parse contigs fasta
    logger.info("Parsing contig fasta file: {}".format(contigs_fasta))
    contig_stats = get_sequence_stats_from_contigs(contigs_fasta)

    # sort and get cumulative length
    contig_stats.fillna(0, inplace=True)
    contig_stats.sort_values(by='Length', ascending=False, inplace=True)
    contig_stats['CumuLength'] = contig_stats.Length.cumsum()
    for col in ['Length', 'ReadCount', 'MaxCov', 'MinCov', 'CumuLength']:
        if col in contig_stats.columns:
            contig_stats[col] = contig_stats[col].astype(int)

    if contig_stats_out:
        logger.info("Writing stats table to: {}".format(contig_stats_out))
        contig_stats.to_csv(contig_stats_out, sep='\t', float_format="%0.2f")
    else:
        return contig_stats

def get_sequence_stats_from_contigs(contigs_fasta):
    """
    Use BioPython parser and GC calculator to get contig lengths and
    GCs from contigs fasta
    """

    # initialize lists
    contigs = []
    lengths = []
    gcs = []

    # loop over fasta records (this is 2-3 times faster than SeqIO.parse)
    # (and only marginally slower than my custom built parser.)
    with open(contigs_fasta, 'r') as contigs_handle:
        for title, sequence in SeqIO.FastaIO.SimpleFastaParser(contigs_handle):
            # parse title with RegEx
            contig = title.split(None, 1)[0]
            length = len(sequence)
            contigs.append(contig)
            lengths.append(length)
            gcs.append(SeqUtils.GC(sequence))

    # convert to DataFrame and return
    return pandas.DataFrame({'Contig': contigs,
                             'Length': lengths,
                             'GC': gcs}).set_index('Contig')


def generate_histogram(contig_stats_file,
                       histogram_file,
                       metric='Length',
                       length_cutoff=2000,
                       log=False,
                       txt_width=75,
                       **kwargs):
    """
    Generate an ASCII histrogram of the selected data
    """
    stats_table = pandas.read_table(contig_stats_file, index_col=0)
    # filter by length, get metric column, replace NA with 0
    values = stats_table \
                .query('Length >= ' + str(length_cutoff)) \
                [metric] \
                .fillna(0)
    with open(histogram_file, 'wt') as hist_handle:
        if values.empty:
            hist_handle.write("There are no contigs longer than "
                              "{length_cutoff}bp".format(**vars()))
        else:
            hist_handle.write(
                """Summary of {metric} for {count} contigs longer than {length_cutoff}bp:
        Min:    {min}
        Max:    {max}
        Mean:   {mean}
        Median: {median}
        StdDev: {std}

    Histogram of {metric} {log}:
    {histogram}
    """.format(metric=metric,
               length_cutoff=length_cutoff,
               count=len(values),
               min=values.min(),
               max=values.max(),
               mean=values.mean(),
               median=values.median(),
               std=values.std(),
               log="(log)" if log else "",
               histogram=ascii_histogram(numpy.histogram(values,
                                                         **kwargs),
                                         log=log,
                                         width=txt_width),
              ))


def get_contig_length_summary_stats(contig_stats, N_levels=[50, 75, 90]):
    """
    uses "length" and "cumulength" columns in contig_stats table
    to quickly get N50 and others for given N_levels (def: 50,75,90)
    """
    if contig_stats.index.empty:
        return {'contig count': 0}

    # tota length is last cumulative length value
    total_length = int(contig_stats['CumuLength'].iloc[-1])

    # init return values with basic data
    N_stats = {'contig count': contig_stats.shape[0],
               'total contig bases': total_length,
              }

    # set up iterator over just these columns
    cumulen_iter = iter(contig_stats[['Length', 'CumuLength']].iterrows())

    # Loop over N's. Since they are sorted, we don't need to restart
    #  the length/cumul_length iterator
    cumulength = 0
    for N in sorted(N_levels):
        # looking for N% of the total length
        target = total_length * N / 100

        # stop when we get there
        while cumulength < target:
            contig, (length, cumulength) = next(cumulen_iter)

        # Save the contig length that got use here
        N_key = "N{0:02d}".format(N)
        N_stats[N_key] = int(length)

    return N_stats
