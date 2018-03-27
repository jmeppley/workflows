###
#
# gene_catalog.py
#
# collection of methods for annotating gene_families
## normalize_coverage(input):
#    process contig and read stats into map from gene to abundance
## 
###
from Bio import SeqIO
import re, os, logging, pandas
from edl import taxon as edltaxon, util, hits as edlhits, blastm8, kegg
from snakemake import logger
from python.common import parse_stats


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

def get_function_trimmed(hit, desc_map):
    """
    Given a RefSeq hit id return a function stripped of boilerplate
    """
    return  func_pref_RE.sub('',
             func_suff_RE.sub('',
              extra_loaction_info_RE.sub(r'\1',
               desc_RE.sub('',desc_map[hit.hit]))))

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
printed_ranks=[edltaxon.ranks[i] for i in [3,7,12,17,21,24]] + ['domain']
major_ranks=[edltaxon.ranks[i] for i in [2,3,7,12,17,21,24,27,28]]
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

class RefSeqGeneAnnotator():
    def __init__(self, db_location,
                 taxid_delim=None,
                 bad_refs=set(),
                 genome_clades={},
                 **kwargs):
        self.rsdb=db_location
        self.parse_db_metadata(taxid_delim)
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
        set_genome_clades( genome_data[column] )

    def parse_db_metadata(self, taxid_delim):
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
        self.hit_translator = edlhits.getHitTranslator(hitStringMap=self.taxid_map, parseStyle=edlhits.ACCS, taxonomy=self.taxonomy, hitsAreObjects=True)

        # descriptions
        rsdb_desc_map = self.rsdb + ".ids"
        self.desc_map = util.parseMapFile(rsdb_desc_map)

    def annotate_genes_rs_prot(self, hit_table, annotation_table):
        with open(annotation_table, 'w') as tsv_out:
            tsv_out.write('Gene\t' \
                           + '\t'.join(printed_ranks) \
                           + '\tfunction\tmin pct ID\thit count' \
                           + '\ttop hit\ttop pct ID\ttop score\ttop desc\n')
            for info in self.generate_gene_annotations_rs_prot(hit_table):
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

    def generate_gene_annotations_rs_prot(self,hit_table):

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
            fns_by_score={}
            for hit in hits:
                f=get_function_trimmed(hit,self.desc_map)
                fns_by_score.setdefault(hit.score,[]).append(f)

            for score in sorted(fns_by_score.keys(),reverse=True):
                # only consider things with useful annotations
                #useful_functions = set(f for f in fns_by_score[score] if uninformative_RE.search(f) is None)
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

            top_hit=hits[0]
            yield (gene,lca_ranked,function,min_pctid,len(hits),top_hit.hit,top_hit.pctid,top_hit.score,self.desc_map.get(top_hit.hit,None))

        logger.info( "Parsed %d hits for %d genes" % (total_hits, total_genes))

ko_hit_org_RE = re.compile(r'^([a-z]+):')
class KeggGeneAnnotator():

    kegg20160201 = '/mnt/lysine/jmeppley/servers/galaxy/py-metagenomics/databases/lastdb/KEGG/KeggGene.pep.20160201/lastdb'

    def __init__(self, db_location=kegg20160201):
        self.keggdb=db_location
        self.parse_db_metadata()
        self.m8_params = edlhits.FilterParams(format=blastm8.BLASTPLUS, top_pct=5., sort='score')
    
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

def normalize_coverages(input):
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
    for cov_file in sorted(input.contig_covs):
        assembly = os.path.dirname(cov_file)
        assembly = re.sub(r'^renamed_dir', '', assembly)
        if assembly != last_assembly:
            if coverages is not None:
                for item in coverages.items():
                    yield item
            coverages = None
            last_assembly = assembly
        _coverages = pandas.read_table(cov_file,
                                       index_col=0,
                                       header=0,
                                       usecols=['Contig', 'MeanCov'],
                                       )['MeanCov'] / read_counts[assembly]
        if coverages is None:
            coverages = _coverages
        else:
            coverages = coverages + _coverages
    for item in coverages.items():
        yield item




