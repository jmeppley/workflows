###
#
# gene_catalog.py
#
# collection of methods for annotating gene_families
## normalize_coverage(input):
#    process contig and read stats into map from gene to abundance
###
import re
import os
import json
import pandas
from Bio import SeqIO
from edl import taxon as edltaxon, util, hits as edlhits, blastm8, kegg
try:
    from python.common import parse_stats
    from snakemake import logger
except:
    # if running as a script
    from common import parse_stats
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

REFSEQ='refseq'
GTDB='GTDB'
KEGG='kegg'

####
# RefSeq
#
# The following methods generate annotations using a REfSeq hit table

# Regular Expressions for parsing functional roles
desc_RE = re.compile(r'\s*\[[^\[]+$')
uninformative_RE=re.compile(r'(:?hypothetical|predicted|conserved\sdomain)\s+protein',)
extra_loaction_info_RE=re.compile(r'^(?:contig)?\d+,\s*\S+_\S+,\s*(\S.+),\s*\d+\s*\-\s*\d+\s*$')
func_pref_RE = re.compile(r'(?:MULTISPECIES|PREDICTED):\s')
func_suff_RE = re.compile(r'\s\((?:plastid|plasmid|chloroplast|chromatophore)\)\s*$')

def get_function_trimmed(hit, desc_map, db_type):
    """
    Given a RefSeq hit id return a function stripped of boilerplate
    """
    if db_type == REFSEQ:
        return  func_pref_RE.sub('',
                 func_suff_RE.sub('',
                  extra_loaction_info_RE.sub(r'\1',
                   desc_RE.sub('',desc_map[hit.hit]))))
    if db_type == GTDB:
        return desc_map[hit.hit].split(None, 1)[-1]
    return "NA"

# drop the accession version suffix
acc_suff_RE = re.compile(r'\.\d+$')
def translate_hit(self, hit):
    hit.hit = acc_suff_RE.sub('',hit.hit)
    return self.hit_translator.translateHit(hit)

# find the least common ancestor of multiple hits
def get_lca(hits, translate=lambda x: [x]):
    orgs = []
    for h in hits:
        orgs.extend(translate(h))
    if len(orgs)==0:
        return None
    lca = orgs.pop()
    while len(orgs)>0:
        lca = lca.getLCA(orgs.pop())
    return lca

# Simplify the list of ranks
printed_ranks=[edltaxon.ranks[i] for i in [3,7,11,17,21,24]] + ['domain']
major_ranks=[edltaxon.ranks[i] for i in [2,3,7,11,17,21,24,27,28]]
def get_major_rank(rank):
    " return the highest major rank below or equal to the given rank "
    rank_index = edltaxon.ranks.index(rank)
    for mr in major_ranks:
        if edltaxon.ranks.index(mr)>=rank_index:
            return mr

def approximate_rank(taxon, use_major_ranks=True):
    ret_rank=None
    if taxon.rank != 'no rank':
        ret_rank = taxon.rank
    else:
        last_rank=edltaxon.ranks[0]
        for rank in edltaxon.ranks:
            if taxon.getAncestorClosestToRank(rank)!=taxon:
                ret_rank = last_rank
                break
            last_rank=rank
        else:
            return "Unknown"

    if use_major_ranks:
        return get_major_rank(ret_rank)
    else:
        return rank

class TaxDBGeneAnnotator():
    def __init__(self, db_location,
                 db_type=REFSEQ,
                 taxid_delim=None,
                 bad_refs=set(),
                 genome_clades={},
                 **kwargs):
        self.rsdb=db_location
        self.parse_db_metadata(taxid_delim, db_type)
        self.set_m8_params(**kwargs)
        self.set_bad_refs(bad_refs)

    def set_m8_params(self, **kwargs):
        kwargs.setdefault('format',blastm8.BLASTPLUS)
        kwargs.setdefault('top_pct',0.)
        self.m8_params = edlhits.FilterParams(**kwargs)

    def set_bad_refs_from_file(self, bad_ref_file, **kwargs):
        """
        Use util.parse_list_to_set() to generate list of bad_refs
        to exclude during annotation
        """
        self.set_bad_refs(util.parse_list_to_set(bad_ref_file, **kwargs))

    def set_bad_refs(self, bad_refs):
        self.bad_refs=bad_refs
        try:
            self.m8_params.bad_refs=bad_refs
        except:
            logger.warn("Can't set bad refs for param")

    def set_genome_clades(self, genome_clades):
        self.genome_clades=genome_clades

    def set_genome_clades_from_file(self,
                                    genome_clade_file,
                                    column='Clade',
                                    **kwargs):
        kwargs.setdefault('index_col',0)
        genome_data=pandas.read_csv(genome_clade_file,**kwargs)
        self.set_genome_clades( genome_data[column] )

    def parse_db_metadata(self, taxid_delim, db_type):
        """
        Given a refseq database in my style, parse the associated map files
        """
        # taxonomy
        rsdb_dir = os.path.split(self.rsdb)[0]
        self.taxonomy = edltaxon.readTaxonomy(rsdb_dir)
        rsdb_taxid_map = self.rsdb + ".tax"
        self.taxid_map = util.parseMapFile(rsdb_taxid_map,
                valueDelim=taxid_delim,
                valueType=int)
        parse_style = edlhits.ACCS if db_type == REFSEQ else edlhits.HITIDS
        self.hit_translator = edlhits.getHitTranslator(hitStringMap=self.taxid_map, parseStyle=parse_style, taxonomy=self.taxonomy, hitsAreObjects=True)

        # descriptions
        rsdb_desc_map = self.rsdb + ".ids"
        self.desc_map = util.parseMapFile(rsdb_desc_map)

    def annotate_genes_rs_prot(self, hit_table, annotation_table,
                               db_type=REFSEQ):
        logger.info("Annotating " + db_type + " taxdb with " + \
                     hit_table + " and " + \
                     annotation_table)
        with open(annotation_table, 'w') as tsv_out:
            tsv_out.write('Gene\t' \
                           + '\t'.join(printed_ranks) \
                           + '\tfunction\tmin pct ID\thit count' \
                           + '\ttop hit\ttop pct ID\ttop score\ttop desc\n')
            for info in self.generate_gene_annotations_rs_prot(hit_table,
                                                               db_type):
                (gene, lca_ranked, function, min_pctid, hit_count, \
                       top_hit, top_pctid, top_score, top_desc) = info
                tsv_out.write("%s\t%s\t%s\t%0.1f\t%d\t%s\t%0.1f\t%0.1f\t%s\n" \
                 % (gene,
                    '\t'.join([lca_ranked.get(r,"") for r in printed_ranks]),
                    function,
                    min_pctid,
                    hit_count,
                    top_hit,
                    top_pctid,
                    top_score,
                    top_desc,
                    ))

    def generate_gene_annotations_rs_prot(self, hit_table, db_type=REFSEQ):
        logger.info("Annotating " + db_type + " taxdb with " + \
                     hit_table)

        species_index = major_ranks.index('species')
        genus_index = major_ranks.index('genus')
        phylum_index = major_ranks.index('phylum')

        total_genes=0
        total_hits=0
        blast_m8 = blastm8.M8Stream(hit_table)
        for gene, hits in blastm8.filterM8Stream(blast_m8,
                self.m8_params,
                returnLines=False):
            hits = list(hits)
            total_hits+=len(hits)
            total_genes+=1
            min_pctid = min([h.pctid for h in hits])
            lca=get_lca(hits, self.hit_translator.translateHit)
            lca_rank = approximate_rank(lca)
            lca_ranked = {'domain': lca.getAncestorClosestToRank('domain').name}
            if lca_rank in major_ranks:
                for r in range(phylum_index,species_index-1,-1):
                    lca_ranked[major_ranks[r]] = lca.getAncestorClosestToRank(major_ranks[r]).name
                    if r<=major_ranks.index(lca_rank):
                        break

            # get a good functional annotation form the best hit(s) (that aren't hypothetical)
            if db_type is not None:
                fns_by_score={}
                for hit in hits:
                    f=get_function_trimmed(hit, self.desc_map, db_type)
                    fns_by_score.setdefault(hit.score,[]).append(f)

                for score in sorted(fns_by_score.keys(),reverse=True):
                    # only consider things with useful annotations
                    useful_functions = set()
                    for f in fns_by_score[score]:
                        if uninformative_RE.search(f) is None:
                            useful_functions.add(f)
                    if len(useful_functions)>0:
                        functions = useful_functions
                        break
                else:
                    functions=['unknown']
                function = ";".join(functions)
            else:
                function = 'NA'


            top_hit=hits[0]
            # description of top hit
            top_desc = self.desc_map[top_hit.hit]
            if db_type == GTDB:
                # GTDB headers are too long, take sp name and func
                top_desc = top_desc.split(";")[-1]
                top_desc = re.sub(r'^s__', '', top_desc)

            yield (gene,lca_ranked,function,min_pctid,len(hits),top_hit.hit,top_hit.pctid,top_hit.score,top_desc)

        logger.info( "Parsed %d hits for %d genes" % (total_hits, total_genes))

ko_hit_org_RE = re.compile(r'^([a-z]+):')
class KeggGeneAnnotator():

    kegg20160201 = '/mnt/lysine/jmeppley/servers/galaxy/py-metagenomics/databases/lastdb/KEGG/KeggGene.pep.20160201/lastdb'

    def __init__(self, db_location=kegg20160201):
        self.keggdb=db_location
        self.parse_db_metadata()
        self.m8_params = edlhits.FilterParams(format=blastm8.BLASTPLUS,
                                              top_pct=5., sort='score')

    def parse_db_metadata(self):
        """
        Given a refseq database in my style, parse the associated map files
        """
        # KO map
        kegg_kos = '%s.kos' % (self.keggdb)
        self.ko_map=kegg.parseLinkFile(kegg_kos)

        # descriptions
        kegg_ids = '%s.ids' % (self.keggdb)
        self.desc_map = util.parseMapFile(kegg_ids)

    def annotate_genes_kg(self, hit_table, annotation_table):
        header='Gene\tKO(s)'
        annot_tuples=self.generate_gene_annotations_kg(hit_table)
        write_tsv(annotation_table, annot_tuples, header=header)

    def generate_gene_annotations_kg(self, hit_table):
        blast_m8 = blastm8.M8Stream(hit_table)
        for read, hits in blastm8.filterM8Stream(blast_m8, self.m8_params, returnLines=False):
            hits = list(hits)
            kos_so_far = []
            org_scores = {}
            for h in hits:
                org=ko_hit_org_RE.search(h.hit).group(1)
                score=h.score
                if org_scores.get(org,-1)>score:
                    # only take top scoring hits from each organism
                    continue
                else:
                    org_scores[org]=score

                kos_so_far.extend(self.ko_map.get(h.hit,[]))
            # collect uniq hits, dropping "None" or None
            uniq_kos=set([k for k in kos_so_far if k not in [None,"None",""]])
            yield (read, ';'.join(uniq_kos))

def write_tsv(out_file, data_tuples, sep='\t', header=None):
    with open(out_file,'w') as out_handle:
        if header is not None:
            out_handle.write(header)
            if not header.endswith('\n'):
                out_handle.write('\n')
        for data_tuple in data_tuples:
            out_handle.write(sep.join(data_tuple)+'\n')

def parse_clusters(cluster_file):
    """
    expects one line per cluster, tab separated:
        cluster_1_rep  member_1_1 member 1_2 ...
        cluster_2_rep  member_2_1 member_2_2 ...
        """
    cluster_dict = {}
    with open(cluster_file) as LINES:
        for line in LINES:
            genes = line.strip().split('\t')
            rep = genes[0]
            for gene in genes:
                cluster_dict[gene] = rep
    return cluster_dict

def parse_bio_clusters(bio_json, out_tab):
    """ convert bio cluster format from vsearch to our syle table """
    with open(bio_json) as BIO:
        data = json.load(BIO)

    # build map from rep to other genes
    clusters = {}
    for row, column, value in data['data']:
        rep = data['rows'][row]['id']
        gene = data['columns'][column]['id']
        if rep == gene:
            # make sure a cluster exists for this rep
            clusters.setdefault(rep, [])
        else:
            # add gene to cluster for this rep
            clusters.setdefault(rep, []).append(gene)

    # write table
    with open(out_tab, 'wt') as TAB:
        for rep, other_genes in clusters.items():
            if len(other_genes) > 0:
                TAB.write("{}\t{}\n".format(rep,
                                            "\t".join(other_genes)))

def parse_mmseq_clusters(mm_tab, out_tab):
    """ reformat mmseqs cluster table to our style 
        input: every line maps rep to member
        output: every non-single cluster listed starting w/rep
    """
    with open(out_tab, 'wt') as TAB:
        with open(mm_tab) as MM:
            prev_rep = None
            gene_count = 0
            cluster_count = 0
            for line in MM:
                rep, gene = line.strip().split('\t')
                if rep != prev_rep:
                    if rep != gene:
                        raise Exception("expected first gene to be same "
                                        " as rep. {} != {}".format(rep, gene))
                    gene_count = 1
                    prev_rep = rep
                    continue
                gene_count += 1
                if gene_count == 2:
                    # start writing, by ending previous line
                    if cluster_count != 0:
                        TAB.write('\n')
                    cluster_count += 1
                    # ... and writing the rep
                    TAB.write(rep)
                # ... add this gene
                TAB.write('\t' + gene)
            # end final cluster
            TAB.write('\n')


def parse_cdhit_clusters(clstr_file, cluster_file):
    """ reformat cdhit's .clstr file into simple table """
    gene_expr = re.compile(r'\s>(\S+)\.\.\.\s\s*(.+)\s*$')
    with open(cluster_file, 'wt') as TAB:
        with open(clstr_file) as CLSTR:
            cluster = []
            cluster_rep = None
            for line in CLSTR:
                if line.startswith('>'):
                    if cluster_rep is not None:
                        TAB.write('{}\t{}\n'.format(
                            cluster_rep,
                            '\t'.join(g for g in cluster)))
                    cluster = []
                    cluster_rep = None
                    continue
                try:
                    gene, alignment = gene_expr.search(line).groups()
                except AttributeError:
                    print("can't parse: \n" + line)
                    raise
                if alignment.strip() == "*":
                    cluster_rep = gene
                else:
                    cluster.append(gene)
            if cluster_rep is not None:
                TAB.write('{}\t{}\n'.format(
                    cluster_rep,
                    '\t'.join(g for g in cluster)))


def merge_cluster_coverages(cluster_file,
                            coverage_tables):
    """
    given a cluster file mapping cluster reps to members
    and a map from assemblies to coverage tables of member genes

    generate a table of cluster coverages by assembly
    """
    cluster_map = parse_clusters(cluster_file)
    cluster_coverages = None
    for assembly, gene_coverage_table in coverage_tables.items():
        gene_coverages = pandas.read_csv(gene_coverage_table,
                                         index_col=0)
        gene_coverages.columns = [assembly,]
        gene_coverages['Cluster'] = [cluster_map.get(g,g) \
                                     for g in gene_coverages.index]
        _cluster_coverages = gene_coverages.groupby('Cluster').agg(sum)
        if cluster_coverages is None:
            cluster_coverages = _cluster_coverages
        else:
            cluster_coverages = \
                cluster_coverages.join(_cluster_coverages,
                                       how='outer')
    return cluster_coverages

# deprecated
def normalize_coverages(input, contig_col='Contig', cov_col='MeanCov'):
    """
    Loop over read stats files and create a normalization factor for each assembly (number of reads/10M)

    input: the input object from snakemake with two member file lists:
    input.read_stats: the cleaned read stats from all assemblies
    input.contig_covs: the contig stats from all assemblies

    Loop over coverage files,
    group by assembly (some assemblies have 2),
    normalize by adjusted number of reads,
    yield (contig, coverage) tuples
    """
    read_counts = {}
    for stats_file in input.read_stats:
        assembly = re.sub(r'/stats/.+$', '', stats_file)
        reads = parse_stats(stats_file)['reads'] / 10000000
        read_counts[assembly] = read_counts.get(assembly, 0) + reads

    last_assembly = None
    coverages = None
    if len(input.contig_covs) == 0:
        raise Exception("Missing contig coverage files!")

    for cov_file in sorted(input.contig_covs):
        assembly = os.path.dirname(cov_file)
        assembly = re.sub(r'^renamed_dir', '', assembly)
        if assembly != last_assembly:
            if coverages is not None:
                for item in coverages.items():
                    yield item
            coverages = None
            last_assembly = assembly
        _coverages = pandas.read_csv(cov_file,
                                     sep='\t',
                                     index_col=0,
                                     header=0,
                                     usecols=[contig_col, cov_col],
                                     )[cov_col] / read_counts[assembly]
        if coverages is None:
            coverages = _coverages
        else:
            coverages = coverages + _coverages
    for item in coverages.items():
        yield item

def main():
    """
    Run as a script:

    gene_catalog.py {hit_table} {lastdb_path} {output_file}
    """
    import argparse
    description = "annotate genes from refseq or kegg hits"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("hit_table", metavar="HIT_TABLE",
                        help="The table of hits in blast format")
    parser.add_argument("lastdb_path", metavar="LASTDB_PATH",
                        help="The path given to lastal")
    parser.add_argument("output_table", metavar="OUT_TABLE",
                        help="The file to write the anntations to")
    parser.add_argument("-t", "--type", default=None, metavar="TYPE",
                        choices=[REFSEQ, KEGG, GTDB],
                        help="The type of database. Either 'refseq'," + \
                        "'GTDB', or 'kegg'")

    arguments = parser.parse_args()

    # try to guess type from DB name/path
    if arguments.type is None:
        if re.search(r'kegg', arguments.lastdb_path,
                     flags=re.I):
            arguments.type = KEGG
        else:
            arguments.type = REFSEQ

    logger.info("Parsing db type: " + arguments.type)

    # annotate!
    if arguments.type == KEGG:
        annotator = KeggGeneAnnotator(arguments.lastdb_path)
        annotator.annotate_genes_kg(arguments.hit_table,
                                    arguments.output_table)
    else:
        annotator = TaxDBGeneAnnotator(arguments.lastdb_path)
        annotator.annotate_genes_rs_prot(arguments.hit_table,
                                         arguments.output_table,
                                         arguments.type)

def process_for_mcl(input_file, fasta_file, output_file,
                    format='last',
                    pctid=.95,
                    minbit=.5):
    """ generates a table of graph edges from an all v all """
    params = blastm8.FilterParams(format=format, pctid=pctid)
    inputm8 = blastm8.M8Stream(input_file)
    # fake all self bits
    self_bits = {r.id: 2*len(r) for r in SeqIO.parse(fasta_file, 'fasta')}
    with open(output_file, 'wt') as OUTPUT:
        for seq, hits in blastm8.filterM8Stream(inputm8, params, returnLines=False):
            for hit in hits:
                if hit.hit == seq:
                    # we've faked the self bits for now
                    continue

                process_hit(hit, OUTPUT, self_bits, minbit)


def process_hit(hit, output_handle, self_bits, minbit):
    """ following mcl parsing used by Anvio, filter hits to feed to mcl """
    bitratio = hit.score / min(self_bits[hit.hit],
                               self_bits[hit.read])
    if bitratio < minbit:
        return
    output_handle.write("{}\t{}\t{}\n".format(hit.hit, hit.read, hit.pctid))


def get_longest_seq(clusters, genes, format='fasta'):
    """
    Given a clsuter file from mcl where each line is a cluster with
    tab separated gene ids.
    And given a fasta file of gene sequences.
    return the ID of the longest gene in each cluster
    """
    gene_lengths = {g.id:len(g) for g in SeqIO.parse(genes, format)}
    with open(clusters) as CLUSTERS:
        for line in CLUSTERS:
            member_genes = line.rstrip().split('\t')
            yield sorted(member_genes, reverse=True,
                         key=lambda g: gene_lengths[g])[0]

if __name__ == '__main__':
    main()
