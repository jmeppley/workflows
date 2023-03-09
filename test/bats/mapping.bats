setup() {
    eval "$(conda shell.bash hook)"
    mkdir -p test/conda/envs
    ENV=mapping
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi

    conda activate $ENV_DIR
}

@test "Map with bwa and a filter" {
    rm -rf test/scratch/map1
    mkdir -p test/scratch/map1
    cd test/scratch/map1
    run bash -c "snakemake -s ../../../scripts/read_mapping.snake --config references_file=../../data/contigs/contigs.aloha.fa sample_glob=../../data/raw_reads/2014_{sample}.fastq -p -k -j 2 > map.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../scripts/read_mapping.snake --config references_file=../../data/contigs/contigs.aloha.fa sample_glob=../../data/raw_reads/2014_{sample}.fastq -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]:0:18}" == "Nothing to be done" ]
}

@test "Map with bwa and simply count " {
    rm -rf test/scratch/map2
    mkdir -p test/scratch/map2
    cd test/scratch/map2
    run bash -c "snakemake -s ../../../scripts/read_mapping.2.snake --config references_file=../../data/contigs/contigs.aloha.fa sample_glob=../../data/raw_reads/2014_{sample}.fastq -p -k -j 2 > map.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../scripts/read_mapping.2.snake --config references_file=../../data/contigs/contigs.aloha.fa sample_glob=../../data/raw_reads/2014_{sample}.fastq -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]:0:18}" == "Nothing to be done" ]
}

@test "Map with kallisto" {
    rm -rf test/scratch/mapk
    mkdir -p test/scratch/mapk
    cd test/scratch/mapk
    run bash -c "snakemake -s ../../../scripts/read_mapping.kallisto.snake --config references_file=../../data/other/assembled_transcripts.fasta sample_glob=../../data/raw_reads/{sample}_merge.fastq -p -k -j 2 > map.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../scripts/read_mapping.kallisto.snake --config references_file=../../data/other/assembled_transcripts.fasta sample_glob=../../data/raw_reads/{sample}_merge.fastq -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]:0:18}" == "Nothing to be done" ]
}
