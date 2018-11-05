"""
Annotation related methods for the assembly workflows:
    drop_all_rna_overlaps: remove genes that overlap [tr]RNA
    filter_annotations: remove annots from "bad" contigs
    filter_and_extract_rRNA: save nonoverlapping rRNA to gff and fna

Other useful functions:
    ranges_intersect: do two start,end pairs overlap by buffer?

"""
import re
from python.common import get_set_from_file
from Bio import SeqIO
from edl.blastm8 import GFF, generate_hits
from snakemake import logger

## RNA and gene filtering

def ranges_intersect(pair1, pair2, buffer=0):
    """ true if ranges entersect by more than buffer """
    return (pair1[0] < pair2[1] - buffer) and (pair1[1] > pair2[0] + buffer)

def gff_location_parser(hit):
    """ get query, start, end from hit """
    return (hit.read, hit.qstart, hit.qend)

def gff_hit_iterator(gff_file, **kwargs):
    """ Loop over gffs individually, not grouped by query(contig) """
    kwargs.setdefault('nonoverlapping', False)
    kwargs.setdefault('sort', 'score')
    for contig, hits in generate_hits(gff_file,
                                      format=GFF,
                                      **kwargs):
        for hit in hits:
            yield hit

def fna_hit_iterator(fasta_file):
    """ loop over fasta records """
    return SeqIO.parse(fasta_file, 'fasta')

# rexp to get contog/start/end from prodigal gene prediction
GENE_REXP = re.compile(
    r"""^(\S+)_\d+\s      # get contig from gene name
        +\#\s+(\d+)\s     # start is first num in #'s
        +\#\s+(\d+)\s     # end is second num in #'s
        +#""", re.X)
# rexp to get contig/start/end from RNA fasta
RNA_GENE_REXP = re.compile(
    r'^(\S+)_[a-z]+RNA_\d+\s.+start=(\d+);end=(\d+);')
def fna_rec_parser(record):
    """ parse start/end from prodigal or RNA fasta header """
    try:
        return GENE_REXP.search(record.description).groups()
    except AttributeError:
        return RNA_GENE_REXP.search(record.description).groups()

def get_mixed_annotations(gff_file):
    """
    Turn the combined annotations gff file into dict mapping
    contigs to list of named gene position tuples
    """
    annots = {}
    for contig, hits in generate_hits(gff_file, format=GFF):
        annot_counts = {t:0 for t in ['CDS', 'rRNA', 'tRNA']}
        contig_annots = annots.setdefault(contig, [])
        for hit in hits:
            gene = get_gene_name(contig, hit, annot_counts)
            start, end = sorted((hit.qstart, hit.qend))
            contig_annots.append((gene, start, end))
    return annots

def get_gene_name(contig, hit, annot_counts):
    """ name gene with count (separate for CDS and [tr]RNA) """
    # is it a CDS or RNA
    if hit.hit_type == 'CDS':
        annot_counts['CDS'] += 1
        suffix = "_{}".format(annot_counts['CDS'])
    elif re.search('rRNA', hit.hit_type) is not None:
        annot_counts['rRNA'] += 1
        suffix = "_rRNA_{}".format(annot_counts['rRNA'])
    elif re.search('tRNA', hit.hit_type) is not None:
        annot_counts['tRNA'] += 1
        suffix = "_tRNA_{}".format(annot_counts['tRNA'])
    else:
        raise Exception("Can't classify GFF line:\n" + hit.line)
    return contig + suffix

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


def drop_all_rna_overlaps(cds_files, rna_files, output_files, buffer):
    """
    Given list of generic annotation files for input
     and list of RNA annotation files for filtering
     and matched list of output annotation files

    print annotations to output files that do not overlap with rna
    """
    # get dict of rRNA locations
    rna_locations = get_annotation_locations(rna_files)

    # filter files
    for input_file in cds_files:
        print("Filtering {}".format(input_file))

        # find output file with same suffix
        output_file = [f for f in output_files \
                            if f.endswith(input_file[-4:])][0]
        drop_rna_overlaps(input_file, output_file, rna_locations,
                          buffer=buffer)


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
    good_contigs = get_set_from_file(good_contig_file)
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
    with open(output_rna_fna, 'wt') as fna_out_handle:
        with open(output_rna_gff, 'wt') as gff_out_handle:
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
                        subsequence = \
                            contig.seq[hit.qstart - 1:hit.qend]
                    else:
                        subsequence = \
                            contig.seq[hit.qend - 1:hit.qstart] \
                            .reverse_complement()
                    if strand == '-':
                        subsequence = subsequence.reverse_complement()

                    fna_out_handle.write(">{seqid}\t{desc}\n{seq}\n"
                                         .format(seqid=gene_name,
                                                 desc=gene_desc,
                                                 seq=str(subsequence)
                                                ))
                    gff_out_handle.write(hit.line)
