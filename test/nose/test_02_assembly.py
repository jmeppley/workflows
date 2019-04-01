" Tests for python/assembly/* "
import os
import python.assembly.annotation as annot
import python.assembly.coverage as cov
import python.assembly.stats as stat
import logging
from edl.blastm8 import Hit, GFF

def test_ranges_intersect():
    pair_1 = (20, 120)
    pair_2 = (90, 200)
    pair_3 = (150, 250)

    assert annot.ranges_intersect(pair_1, pair_2)
    assert annot.ranges_intersect(pair_2, pair_3)
    assert not annot.ranges_intersect(pair_1, pair_3)
    assert annot.ranges_intersect(pair_1, pair_2, 20)
    assert annot.ranges_intersect(pair_1, pair_2, buffer=20)
    assert not annot.ranges_intersect(pair_1, pair_2, 40)
    assert not annot.ranges_intersect(pair_1, pair_2, buffer=40)
    assert annot.ranges_intersect(pair_2, pair_3, 40)
    assert not annot.ranges_intersect(pair_2, pair_3, 60)

SPADES_GFF_LINE = \
    "spades_c1	Prodigal_v2.6.3	CDS	721	822	1.6	-	0	ID=1_1;partial=00;start_type=ATG;rbs_motif=AGxAGG/AGGxGG;rbs_spacer=5-10bp;gc_cont=0.422;conf=59.08;score=1.60;cscore=-2.64;sscore=4.24;rscore=3.31;uscore=0.36;tscore=1.07;\n"
RNA_GFF_LINE = \
    "spades_c104	cmsearch	LSU_rRNA_bacteria	2	278	135.6	+	.	Target=RF02541 1162 1453\n"

def test_gff_parser():
    gff_hit = Hit(SPADES_GFF_LINE, GFF)
    assert (annot.gff_location_parser(gff_hit) ==
            ('spades_c1', 721, 822))
    assert (annot.gff_location_parser(Hit(RNA_GFF_LINE, GFF)) ==
            ('spades_c104', 2, 278))


SPADES_FFN_DESC = \
    "spades_c1_1 # 721 # 822 # -1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=AGxAGG/AGGxGG;rbs_spacer=5-10bp;gc_cont=0.422"
RNA_FNA_DESC = \
    "spades_c1_rRNA_1	source=cmsearch;type=LSU_rRNA_bacteria;score=1635.2;strand=-;start=2208;end=2;Target=RF02541 369 2558contig_length=2209;contig_cov=10.605668"

def test_fna_parser():
    assert (annot.GENE_REXP.search(SPADES_FFN_DESC).groups() ==
            ("spades_c1", "721", "822"))
    assert (annot.RNA_GENE_REXP.search(RNA_FNA_DESC).groups() ==
            ("spades_c1", "2208", "2"))

def test_gff_functions():
    rna_gff_file = "test/bats/outputs/nose/contigs.all.annotations.gff"
    hit_locations = annot.get_annotation_locations([rna_gff_file,])
    assert isinstance(hit_locations, dict)
    contig, locations = next(iter(hit_locations.items()))
    assert isinstance(contig, str)
    logging.debug("%s: %r", contig, locations)
    location = next(iter(locations))
    assert len(location) == 2
    assert isinstance(location[0], int)
    assert isinstance(location[1], int)

def test_get_gene_name():
    gff_hit = Hit(SPADES_GFF_LINE, GFF)
    assert (annot.get_gene_name(gff_hit, {'CDS': 2}) ==
            'spades_c1_3')
    gff_hit = Hit(RNA_GFF_LINE, GFF)
    assert (annot.get_gene_name(gff_hit, {'rRNA': 1}) ==
            'spades_c104_rRNA_2')


## coverage
DEPTHS_LINES = ["contig_c1\t5\t1\n",
                "contig_c1\t6\t1\n",
                "contig_c1\t7\t2\n",
                "contig_c1\t8\t2\n",
                "contig_c1\t9\t2\n",
                "contig_c1\t10\t2\n",
                "contig_c2\t4\t1\n",
                "contig_c2\t5\t1\n",
                "contig_c2\t6\t1\n"]

CONTIG_LENGTHS = {'contig_c1':20, 'contig_c2':9, 'contig_c3':5}

def test_contig_coverage_table():
    g = cov.contig_depths_generator(DEPTHS_LINES)
    contig, depths = next(g)
    assert contig == 'contig_c1'
    assert len(depths) == 6
    contig, depths = next(g)
    assert len(depths) == 3
    coverage = list(cov._insert_zeros(depths, 9))
    assert len(coverage) == 9
    assert sum(coverage) == 3

    coverages = cov.get_contig_coverage_table(DEPTHS_LINES, None)
    logging.debug(str(coverages))
    assert coverages.shape == (2, 6)
    assert coverages.loc['contig_c1', 'MeanCov'] == 1.0
    assert coverages.loc['contig_c2', 'Q2Q3Cov'] == 0.5

    coverages = cov.get_contig_coverage_table(DEPTHS_LINES,
                                              CONTIG_LENGTHS)
    logging.debug(str(coverages))
    assert coverages.shape == (3, 6)
    assert coverages.loc['contig_c2', 'Q2Q3Cov'] == 0.2

def test_contig_stats():
    contigs_file = 'test/data/contigs/contigs.100.fa'
    stats_table = stat.get_sequence_stats_from_contigs(contigs_file)
    logging.debug(str(stats_table))
    assert stats_table.shape == (100, 2)
    assert '{:0.3f}'.format(stats_table.loc['k99_227', 'GC']) == '36.908'
    assert stats_table.loc['k99_227', 'Length'] == 2723
    assert stats_table.loc['k99_7348', 'Length'] == 9664
    assert stats_table.loc['k99_10078', 'Length'] == 10609

    stats_table_2 = stat.get_contig_stats(contigs_file)
    assert stats_table_2.shape == (100, 3)
    assert '{:0.3f}'.format(stats_table_2.loc['k99_227', 'GC']) == '36.908'
    assert stats_table_2.loc['k99_227', 'Length'] == 2723
    assert stats_table_2.loc['k99_7348', 'Length'] == 9664
    assert stats_table_2.loc['k99_10078', 'Length'] == 10609
    last_contig = stats_table_2.index[-1]
    first_contig = stats_table_2.index[0]
    assert stats_table_2.loc[first_contig, 'Length'] == \
            stats_table_2.loc[first_contig, 'CumuLength']
    assert stats_table_2.loc[last_contig, 'CumuLength'] == \
            stats_table_2.Length.sum()


    summary_stats = stat.get_contig_length_summary_stats(stats_table_2,
                                                         N_levels=[50, 80])

    assert summary_stats['contig count'] == 100
    assert summary_stats['total contig bases'] == stats_table.Length.sum()
    logging.debug(summary_stats)
    assert summary_stats['N80'] == 2989

