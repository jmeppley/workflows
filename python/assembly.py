from Bio import SeqIO, SeqUtils
from edl.util import ascii_histogram
import numpy
import pandas
from snakemake import logger


def get_contig_stats(contigs_fasta,
                     contig_stats_out):
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

    logger.info("Writing stats table to: {}".format(contig_stats_out))
    contig_stats.to_csv(contig_stats_out, sep='\t', float_format="%0.2f")


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
    values = stats_table.query('Length >= ' + str(length_cutoff))[metric]
    with open(histogram_file, 'wt') as HOUT:
        HOUT.write("""Summary of {metric} for {count} contigs longer than {length_cutoff}bp:
    Min:    {min}
    Max:    {max}
    Mean:   {mean}
    Median: {median}
    StdDev: {std}

Histogram of {metric} {log}:
{histogram}
"""              .format(
                        metric=metric,
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
                   )
                  )


def get_coverage_stats(contig_depth_file,
                       contig_read_counts_file,
                       contig_stats_out):
    """
    Generate table of read counts and mapped coverage from samtools output
    """
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
    mapping_depth_table = get_samtool_depth_table(contig_depth_file)
    contig_stats = mapping_depth_table.join(read_count_table, how='left').fillna(0)

    for col in ['Length', 'ReadCount', 'MaxCov', 'MinCov', 'CumuLength']:
        if col in contig_stats.columns:
            contig_stats[col] = contig_stats[col].astype(int)

    logger.info("Writing coverage table to: {}".format(contig_stats_out))
    contig_stats.to_csv(contig_stats_out, sep='\t', float_format="%0.2f")


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
    with open(contigs_fasta, 'r') as CF:
        for title, sequence in SeqIO.FastaIO.SimpleFastaParser(CF):
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

def get_samtool_depth_table(depth_file):
    """
    Calculate coverage stats for each contig in an assembly

    Params:
     depth_file: output file from the command:
                    `samtools depth reads.v.contigs.bam`
                 this is a 3 column file with one line per base.
                 columns are:
                     'contig_id base_index base_depth'

    Returns:
     pandas.DataFrame with one row per contig and the three following columns:
            Contig, MeanCov, MedCov, MaxCov, MinCov, StdCov
            """
    with open(depth_file, 'r') as DEPTHS:
        return get_samtool_depth_table_from_handle(DEPTHS)


def get_samtool_depth_table_from_handle(depth_stream):
    """
    Calculate coverage stats for each contig in an assembly

    Params:
     depth_stream: output file from the command:
                    `samtools depth reads.v.contigs.bam`

                    passed as an open file-like object (aka a file handle)
                 this is a 3 column file with one line per base.
                 columns are:
                     'contig_id base_index base_depth'

    Returns:
     pandas.DataFrame with one row per contig and the three following columns:
            Contig, MeanCov, MedCov, MaxCov, MinCov, StdCov
            """

    # reading into lists is a fast way to build a big DataFrame
    contigs, av_covs, mn_covs, mx_covs, md_covs, std_covs = \
            [], [], [], [], [], []

    # loop over contig bases
    current_contig = None
    depths = []
    for line in depth_stream:
        contig, base, depth = line.split()
        depth = int(depth)
        if contig != current_contig:
            if current_contig is not None:
                # end of contig, save numbers
                contigs.append(current_contig)
                depth_array = numpy.array(depths)
                av_covs.append(depth_array.mean())
                mn_covs.append(depth_array.min())
                mx_covs.append(depth_array.max())
                md_covs.append(numpy.median(depth_array))
                std_covs.append(depth_array.std())
            depths = []
            current_contig = contig

        # update contig numbers with current base
        depths.append(depth)

    # end of final contig, save numbers
    contigs.append(current_contig)
    depth_array = numpy.array(depths)
    av_covs.append(depth_array.mean())
    mn_covs.append(depth_array.min())
    mx_covs.append(depth_array.max())
    md_covs.append(numpy.median(depth_array))
    std_covs.append(depth_array.std())

    #return pandas.DataFrame(
    #    [contigs, av_covs, mx_covs, mn_covs, std_covs, md_covs],
    #    columns=['Contig', 'MeanCov', 'MaxCov', 'MinCov', 'StdCov', 'MedCov'],
    #                       ).set_index('Contig')
    table = pandas.DataFrame([contigs,
                              av_covs,
                              mx_covs,
                              mn_covs,
                              std_covs,
                              md_covs]).T
    table.columns=['Contig',
                   'MeanCov',
                   'MaxCov',
                   'MinCov',
                   'StdCov',
                   'MedCov']
    return table.set_index('Contig')


def get_contig_length_summary_stats(contig_stats, N_levels=[50, 75, 90]):
    """
    uses "length" and "cumulength" columns in contig_stats table
    to quickly get N50 and others for given N_levels (def: 50,75,90)
    """
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
    for N in sorted(N_levels):
        # looking for N% of the total length
        target = total_length * N / 100

        # stop when we get there
        cumulength = 0
        while cumulength < target:
            contig, (length, cumulength) = next(cumulen_iter)

        # Save the contig length that got use here
        N_key = "N{0:02d}".format(N)
        N_stats[N_key] = int(length)

    return N_stats



