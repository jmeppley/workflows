setup() {
    ENV=illumina.qc
    ENV_DIR=`pwd`/test/conda/envs/illumina.qc
    ENV_FILE=test/conda/illumina.qc.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet
    fi
    source activate $ENV_DIR

    rm -rf test/scratch/qc
    mkdir -p test/scratch/qc
    cd test/scratch/qc
    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R1_001.head40k.fastq aloha1b.R1.fastq
    ln -s ../../data/raw_reads/2014_ALOHA_XVII_1-1B_S1_R2_001.head40k.fastq aloha1b.R2.fastq
}

@test "Prep a pair of scriptseq files with pear" {
    run bash -c "snakemake -s ../../../test.snake -p --config workflow=qc/pear.snake --verbose aloha1b.scripseq.ATCACG.trim_adapt.joined.fastq > pear.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Prep a pair of scriptseq files with flash" {
    run bash -c "snakemake -s ../../../test.snake -p --config workflow=qc/flash.snake --verbose aloha1b.scripseq.ATCACG.trim_adapt.joined.fastq > flash.join.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Prep a pair of scriptseq files with pandaseq" {
    run bash -c "snakemake -s ../../../test.snake -p --config workflow=qc/pandaseq.snake --verbose aloha1b.scripseq.ATCACG.trim_adapt.joined.fastq > pandaseq.join.log 2>&1"
    [ "$status" -eq 0 ]
}

