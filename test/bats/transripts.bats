setup() {
    ENV=salmon
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    SALMONENV=$ENV_DIR
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf "$ENV_DIR"
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi

    ENV=transcripts
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi
    # suppress warnings
    source activate $ENV_DIR 2>/dev/null || true
    # quit if activation failed
    if [ "$CONDA_PREFIX" != "$ENV_DIR" ]; then
        echo "ERROR activating cond env $ENV_DIR"
        exit 1
    fi
}

@test "assemble single transcriptome with spades" {
    rm -rf test/scratch/transcripts
    mkdir -p test/scratch/transcripts
    cd test/scratch/transcripts

    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config salmon_env=$SALMONENV --verbose > assembly.log 2>&1"
}
