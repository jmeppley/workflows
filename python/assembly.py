"""
Functions used in the assembly workflows
"""
import re
from collections import defaultdict
import numpy
import pandas
from Bio import SeqIO, SeqUtils
from edl.util import ascii_histogram
from edl.blastm8 import GFF, generate_hits
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

def get_samtool_depth_table(depth_file, fasta_file):
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
                      SeqIO.parse(fasta_file, 'fasta')}

    with open(depth_file) as depth_handle:
        return DepthParser(depth_handle).get_data_frame(contig_lengths)

class DepthParser():
    """
    Calculate coverage stats for each contig in an assembly

    Requires open file-like handle to `samtools depths` output
                 this is a 3 column file with one line per base.
                 columns are:
                     'contig_id base_index base_depth'
    contg_lengths: (optional) dict from contig to length
                This is un-needed if depths includes bases with 0 depth

    """

    def __init__(self, depth_file_handle):
        """ initialze parser with open file handle on depths """
        self.handle = depth_file_handle

    # always look ahead to see if contig changes
    next_entry = None

    def __iter__(self):
        """ return an iterator over contigs"""
        try:
            self.next_entry = _parse_depth_line(next(self.handle))
        except StopIteration as e:
            pass
        return self

    def __next__(self):
        """ if there is a next line, return an iterator over depths """
        if self.next_entry is not None:
            return self.next_entry[0], self._generate_depths_for_next_contig()
        else:
            raise StopIteration

    def _generate_depths_for_next_contig(self):
        """ yield lines until the contig changes """
        # remember which contig we're on
        this_contig = self.next_entry[0]
        yield self.next_entry[1]

        # yield lines until we see a new contig
        for line in self.handle:
            self.next_entry = _parse_depth_line(line)
            contig = self.next_entry[0]
            if contig == this_contig:
                yield self.next_entry[1]
            else:
                break
        else:
            # we reached the end
            self.next_entry = None

    columns = ['Contig',
               'MedCov',
               'Q2Q3Cov',
               'MeanCov',
               'StdCov',
               'MinCov',
               'MaxCov',
              ]

    def get_data_frame(self, contig_lengths=None):
        """
        Uses this parser to generate a DataFrame of coverage stats by contig
        Pass in dict of contig legths if 0 depths are excluded from depth file
        """
        if contig_lengths is None:
            contig_lengths = defaultdict(lambda: 0)
        return pandas.DataFrame(self._contig_data_generator(contig_lengths),
                                columns=self.columns) \
                .set_index('Contig')

    def _contig_data_generator(self, contig_lengths):
        for contig, depths in self:
            coverage = numpy.array(list(_insert_zeros(depths, contig_lengths[contig])))
            yield (contig,
                   numpy.median(coverage),
                   q2q3_mean(coverage),
                   coverage.mean(),
                   coverage.std(),
                   coverage.min(),
                   coverage.max())


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
    q = int(len(coverage) / 4)
    return coverage[q:-q].mean()

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

## RNA and gene filtering

def ranges_intersect(pair1, pair2, buffer=0):
    return (pair1[0] < pair2[1] - buffer) and (pair1[1] > pair2[0] + buffer)

def gff_location_parser(hit):
    return (hit.read, hit.qstart, hit.qend)

def gff_hit_iterator(gff_file, nonoverlapping=False):
    for c, h in generate_hits(gff_file,
                              format=GFF,
                              sort='score',
                              nonoverlapping=nonoverlapping):
        for hit in h:
            yield hit

def fna_hit_iterator(fasta_file):
    return SeqIO.parse(fasta_file, 'fasta')

gene_rexp = re.compile(r'^>?((?<!>)\S+)_\d+\s+#\s+(\d+)\s+#\s+(\d+)\s+#')
rna_gene_rexp = \
        re.compile(r'^>?((?<!>)\S+)_\d+\s.+start=(\d+);end=(\d+);')
def fna_rec_parser(record):
    """ parse start/end from prodigal fasta header """
    try:
        return gene_rexp.search(record.description).groups()
    except AttributeError:
        return rna_gene_rexp.search(record.description).groups()

def get_annotation_locations(gff_files):
    """ Given some gff files, build dict of hit locations """
    if isinstance(gff_files, str):
        gff_files = [gff_files, ]

    hit_locations = {}
    for gff_file in gff_files:
        for gene in gff_hit_iterator(gff_file, nonoverlapping=True):
            contig, start, end = gff_location_parser(gene)
            new_rna_pos = sorted([int(p) for p in (start, end)])
            hit_locations.setdefault(contig, []).append(new_rna_pos)

    return hit_locations

def get_good_contig_list(good_contigs):
    """ parse contig list into set """
    with open(good_contigs) as contig_handle:
        good_contigs = set(c.strip() for c in contig_handle.readlines())
    return good_contigs


def get_file_handlers(input_file):
    """ based on extension return tools for looping over annots """
    # set up parsing by file type
    file_type = input_file[-3:]
    if file_type.startswith('f'):
        iterator = fna_hit_iterator
        record_parser = fna_rec_parser
        formatter = lambda g: g.format('fasta')
    else:
        iterator = gff_hit_iterator
        record_parser = gff_location_parser
        formatter = lambda l: l.line
    return iterator, record_parser, formatter


def filter_annotations(good_contig_file, input_file, output_file):
    """ Filter the annotations in the given file (input_file) using
    the supplied list of good contigs """
    good_contigs = get_good_contig_list(good_contig_file)
    iterator, record_parser, formatter = get_file_handlers(input_file)

    # Only keep genes in good contigs
    in_genes = 0
    out_genes = 0
    with open(output_file, 'wt') as out_handle:
        for gene in iterator(input_file):
            in_genes += 1
            contig, start, end = record_parser(gene)
            if contig in good_contigs:
                out_genes += 1
                out_handle.write(formatter(gene))
    logger.debug("Kept {} of {} genes from {}" \
                 .format(out_genes, in_genes, input_file))


def drop_rna_overlaps(input_file, output_file,
                      rna_locations, buffer=0):
    """
    For gene in input file (fasta or gff), print out genes that
    do not overlap rna_locations to output_file
    """
    iterator, record_parser, formatter = get_file_handlers(input_file)

    # Only keep genes that don't overlap
    in_genes = 0
    out_genes = 0
    with open(output_file, 'wt') as out_handle:
        for gene in iterator(input_file):
            in_genes += 1
            contig, start, end = record_parser(gene)
            gene_pos = sorted([int(p) for p in (start, end)])
            for rna_pos in rna_locations.get(contig, []):
                if ranges_intersect(gene_pos, rna_pos, buffer):
                    break
            else:
                out_genes += 1
                out_handle.write(formatter(gene))
    logger.debug("Kept {} of {} genes from {}" \
                 .format(out_genes, in_genes, input_file))


def filter_and_extract_rRNA(raw_gff, contigs_fasta,
                            output_rna_fna,
                            output_rna_gff,
                            molecule,
                            buffer=0):
    """
    Take raw GFF of rRNA hits, and output GFF and fna of non-overlapping
    """
    filter_args = {'sort': 'score',
                   'nonoverlapping': buffer}
    gff_hits = {c: list(h) for c, h in generate_hits(raw_gff,
                                                     format=GFF,
                                                     **filter_args)}
    with open(output_rna_fna, 'wt') as FNA_OUT:
        with open(output_rna_gff, 'wt') as GFF_OUT:
            # loop over contigs
            for contig in SeqIO.parse(contigs_fasta, 'fasta'):
                # get naming informatoion from gff line and contig
                # from spades:
                m = re.search(r'length_(\d+)_cov_([0-9.]+)',
                              contig.description)
                if m:
                    length, coverage = m.groups()
                    contig_info_string = \
                            "contig_length={};contig_cov={}" \
                                        .format(length, coverage)
                else:
                    # megahit
                    m = re.search(r'\blen=(\d+)\b', contig.description)
                    if m:
                        contig_info_string = "contig_length=" + m.group(1)
                    else:
                        contig_info_string = ""

                # loop over hits
                for i, hit in enumerate(gff_hits.get(contig.id, [])):
                    source, feature_type, start, end, score, strand = \
                        hit.line.split('\t')[1:7]

                    # name gene with contig name and index
                    gene_name = contig.id + "_{mol}_{n}" \
                                                .format(n=i+1,
                                                        mol=molecule)
                    # put everything else in the description
                    gene_desc =\
                        ("source={source};type={feature_type};"
                         "score={score};"
                         "strand={strand};start={start};end={end};") \
                        .format(source=source, start=start, end=end,
                                feature_type=feature_type, score=score,
                                strand=strand) + \
                        hit.hitDesc + \
                        contig_info_string

                    # reverse comp seq if on other strand
                    if hit.qstart <= hit.qend:
                        subsequence = contig.seq[hit.qstart - 1:hit.qend]
                    else:
                        subsequence = contig.seq[
                            hit.qend - 1:hit.qstart].reverse_complement()
                    if strand == '-':
                        subsequence = subsequence.reverse_complement()

                    FNA_OUT.write(">{seqid}\t{desc}\n{seq}\n"
                                            .format(seqid=gene_name,
                                                    desc=gene_desc,
                                                    seq=str(subsequence)
                                                    ))
                    GFF_OUT.write(hit.line)


