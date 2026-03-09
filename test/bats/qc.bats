setup() {
    eval "$(conda shell.bash hook)"
    mkdir -p test/conda/envs
    ENV=illumina.qc
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --quiet > test/conda/envs/.create.$ENV 2>&1
    fi
    conda activate $ENV_DIR
}

@test "Prep a pair of scriptseq files with pear" {
    rm -rf test/scratch/qc.pear
    mkdir -p test/scratch/qc.pear
    cd test/scratch/qc.pear

    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R1_001.head40k.fastq aloha1b.R1.fastq
    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R2_001.head40k.fastq aloha1b.R2.fastq
    run bash -c "snakemake -j 10 -s ../../../meta.snake -p --config workflow=qc/pear.snake --verbose aloha1b.scripseq.ATCACG.trim_adapt.joined.fastq > pear.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Prep two pair of scriptseq files with flash" {
    rm -rf test/scratch/qc.flash
    mkdir -p test/scratch/qc.flash
    cd test/scratch/qc.flash

    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config joining_program=flash discover_fastx_for_stats=True > flash.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Split reads into rRNA or not" {
    rm -rf test/scratch/qc.rna
    mkdir -p test/scratch/qc.rna
    cd test/scratch/qc.rna

    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R2_001.head40k.fastq aloha1b.R2.fastq
    run bash -c "snakemake -j 10 -s ../../../meta.snake -p --config workflow=qc/sort.rna.snake file_root=aloha1b.R2 -k sort_rna_default_all > sort.rna.log 2>&1"
    [ "$status" -eq 0 ]
}
