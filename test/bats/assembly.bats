setup() {
    eval "$(conda shell.bash hook)"
    mkdir -p test/conda/envs
    ENV=assembly
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi
    # suppress warnings
    conda activate $ENV_DIR 2>/dev/null || true
    # quit if activation failed
    if [ "$CONDA_PREFIX" != "$ENV_DIR" ]; then
        echo "ERROR activating cond env $ENV_DIR"
        exit 1
    fi
}

@test "assemble single sample with spades" {
    rm -rf test/scratch/spades
    mkdir -p test/scratch/spades
    cd test/scratch/spades

    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --config discover_fastx_for_stats=True --configfile ../../data/configs/spades.yaml --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    [ -e stats/fourk.renamed.interleaved.noadapt.nophix.fastq.stats ]
    [ -e stats/fourk.renamed.interleaved.noadapt.fastq.hist ]
    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/spades.yaml -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
    run bash -c "head -1 fortyk.renamed.R1.tsv | perl -pe 's/\s//'g"
    [ "${lines[0]}" == "NS500472:37:HMMVTBGXX:1:11101:14387:1050fortyk_1" ]
}

@test "assemble clean reads with spades" {
    rm -rf test/scratch/spadesc
    mkdir -p test/scratch/spadesc
    cd test/scratch/spadesc

    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --config discover_fastx_for_stats=True --configfile ../../data/configs/spades.clean.yaml --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    [ -e stats/contigs.filter.fasta.hist ]
    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/spades.clean.yaml -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "assemble clean reads with spades applying bfc" {
    rm -rf test/scratch/spadesec
    mkdir -p test/scratch/spadesec
    cd test/scratch/spadesec

    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --config discover_fastx_for_stats=True --configfile ../../data/configs/spades.econly.yaml --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    [ -e stats/contigs.filter.fasta.hist ]
    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/spades.econly.yaml -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "assemble filtered reads with spades" {
    rm -rf test/scratch/filter
    mkdir -p test/scratch/filter
    cd test/scratch/filter

    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --config discover_fastx_for_stats=True --configfile ../../data/configs/spades_filtered.yaml --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    [ -e stats/fortyk.renamed.interleaved.noadapt.nophix.fastq.stats ]
    [ -e stats/fortyk.renamed.interleaved.noadapt.fastq.hist ]
    [ -e stats/contigs.filter.fasta.stats ]
    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/spades_filtered.yaml -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "assemble single sample with megahit" {
    rm -rf test/scratch/megahit
    mkdir -p test/scratch/megahit
    cd test/scratch/megahit

    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/megahit.yaml --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/megahit.yaml -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Assemble with megahit and fastp" {
    rm -rf test/scratch/fp-megahit
    mkdir -p test/scratch/fp-megahit
    cd test/scratch/fp-megahit

    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/megahit.yaml --config cleaning_protocol=assembly_fastp --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/megahit.yaml --config cleaning_protocol=assembly_fastp -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Assemble with megahit and iutils" {
    rm -rf test/scratch/iu-megahit
    mkdir -p test/scratch/iu-megahit
    cd test/scratch/iu-megahit

    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/megahit.yaml --config cleaning_protocol=anvio --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/megahit.yaml --config cleaning_protocol=anvio -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Prep two pair of files for assembly" {
    rm -rf test/scratch/qc.assembly
    mkdir -p test/scratch/qc.assembly
    cd test/scratch/qc.assembly

    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config cleaning_protocol=assembly remove_rna=False discover_fastx_for_stats=True > assembly.clean.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --config cleaning_protocol=assembly remove_rna=False discover_fastx_for_stats=True -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}
