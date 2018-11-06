" Tests for python/assembly/* "
import os
import python.assembly.annotation as annot
import python.assembly.coverage as cov
import python.assembly.stats as stat
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
            ('spades_c1', 712, 822))
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
    assert isinstance(next(iter(hit_locations)), str)
    assert isinstance(next(iter(hit_locations.values())), tuple)

def test_get_gene_name():
    gff_hit = Hit(SPADES_GFF_LINE, GFF)
    assert (annot.get_gene_name(gff_hit, {'CDS', 2}) ==
            'spades_c1_3')
    gff_hit = Hit(RNA_GFF_LINE, GFF)
    assert (annot.get_gene_name(gff_hit, {'rRNA', 1}) ==
            'spades_c104_rRNA_2')
