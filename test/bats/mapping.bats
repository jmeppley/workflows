setup() {
    if [ "$?" -eq "1" ]; then
        unset conda
        export PATH=$_CONDA_ROOT:$PATH
    fi

    mkdir -p test/conda/envs
    ENV=mapping
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi

    source activate $ENV_DIR
}

@test "Map with bwa and a filter" {
    rm -rf test/scratch/map1
    mkdir -p test/scratch/map1
    cd test/scratch/map1
    run bash -c "snakemake -s ../../../scripts/read_mapping.snake --config references_file=../../data/other/all_genes.clustered.faa sample_glob=../../data/raw_reads/2014_{sample}.fastq -p -k -j 2 > map.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../scripts/read_mapping.snake --config references_file=../../data/other/all_genes.clustered.faa sample_glob=../../data/raw_reads/2014_{sample}.fastq -p -k -j 2 -n 2>&1 | grep Nothing"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}
