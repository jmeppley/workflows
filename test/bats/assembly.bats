setup() {
    ENV=assembly
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

@test "assemble single sample with spades" {
    rm -rf test/scratch/spades
    mkdir -p test/scratch/spades
    cd test/scratch/spades

    run bash -c "snakemake -j 10 -s ../../../assembly.metagenomic.snake -p --configfile ../../data/configs/spades.yaml --verbose > assembly.log 2>&1"
}

@test "Prep two pair of files for assembly" {
    rm -rf test/scratch/qc.assembly
    mkdir -p test/scratch/qc.assembly
    cd test/scratch/qc.assembly

    run bash -c "snakemake -j 10 -s ../../../clean.illumina.snake -p --configfile ../../data/configs/illumina.yaml --verbose --config cleaning_protocol=assembly remove_rna=False discover_fastx_for_stats=True > assembly.clean.log 2>&1"
    [ "$status" -eq 0 ]
}

