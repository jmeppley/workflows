import re
from Bio import SeqIO
from edl.blastm8 import GFF, generate_hits
from snakemake import logger

## RNA and gene filtering

def ranges_intersect(pair1, pair2, buffer=0):
    return (pair1[0] < pair2[1] - buffer) and (pair1[1] > pair2[0] + buffer)

def gff_location_parser(hit):
    return (hit.read, hit.qstart, hit.qend)

def gff_hit_iterator(gff_file, **kwargs):
    kwargs.setdefault('nonoverlapping', False)
    kwargs.setdefault('sort', 'score')
    for c, h in generate_hits(gff_file,
                              format=GFF,
                              **kwargs):
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

def get_mixed_annotations(annot_file):
    """
    Turn the combined annotations gff file into dict mapping
    contigs to list of named gene position tuples
    """
    annots = {}
    for contig, hits in generate_hits(gff_file, format='gff'):
        annot_counts = {t:0 for t in ['CDS','rRNA','tRNA']}
        contig_annots = annots.setdefault(contig, [])
        for hit in hits:
            gene = get_gene_name(contig, hit, annot_counts)
            start, end = sorted(hit.qstart, h.qend)
            contig_annots.append((gene, start, end))
    return annots

def get_gene_name(contig, hit, annot_counts):
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


