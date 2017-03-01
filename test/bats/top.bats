setup() {
    ENV=illumina.qc
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet
    fi
    source activate $ENV_DIR
}

@test "Find top hits in a database" {
    rm -rf test/scratch/top
    mkdir -p test/scratch/top
    cd test/scratch/top
    run bash -c "snakemake -s ../../../annotation.tophits.snake --configfile ../../data/configs/reads.top.RS.yaml -p -k -j 20 > top.annot.log 2>&1"
    [ "$status" -eq 0 ]
}

@test "Find top hits in a database with QC" {
    rm -rf test/scratch/topqc
    mkdir -p test/scratch/topqc
    cd test/scratch/topqc
    run bash -c "snakemake -s ../../../annotation.tophits.snake --configfile ../../data/configs/reads.top.QC.RS.yaml -p -k -j 20 > top.annot.log 2>&1"
    [ "$status" -eq 0 ]
}

