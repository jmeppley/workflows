setup() {
    eval "$(conda shell.bash hook)"
    mkdir -p test/conda/envs
    ENV=transcripts
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
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

@test "assemble single transcriptome with spades, map with kallisto" {
    rm -rf test/scratch/transcripts.kallisto
    mkdir -p test/scratch/transcripts.kallisto
    cd test/scratch/transcripts.kallisto

    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config long_ssu_length=150 long_lsu_length=150 --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config long_ssu_length=150 long_lsu_length=150 -n 2>&1 | grep '^Nothing to be done'"
    [ "$status" -eq 0 ]
    [ "${lines[0]:0:18}" == "Nothing to be done" ]
}

@test "assemble single transcriptome with spades, map with salmon" {
    rm -rf test/scratch/transcripts.salmon
    mkdir -p test/scratch/transcripts.salmon
    cd test/scratch/transcripts.salmon

    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config long_ssu_length=150 long_lsu_length=150 transcript_mapper=salmon --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config long_ssu_length=150 long_lsu_length=150 transcript_mapper=salmon -n 2>&1 | grep '^Nothing to be done'"
    [ "$status" -eq 0 ]
    [ "${lines[0]:0:18}" == "Nothing to be done" ]
}

@test "assemble single transcriptome with spades, map with bwa" {
    rm -rf test/scratch/transcripts.bwa
    mkdir -p test/scratch/transcripts.bwa
    cd test/scratch/transcripts.bwa

    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config long_ssu_length=150 long_lsu_length=150 map_clean_reads=True transcript_mapper=bwa --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config long_ssu_length=150 long_lsu_length=150 map_clean_reads=True transcript_mapper=bwa -n 2>&1 | grep '^Nothing to be done'"
    [ "$status" -eq 0 ]
    [ "${lines[0]:0:18}" == "Nothing to be done" ]
}

@test "assemble single transcriptome with spades, don't use bfc" {
    rm -rf test/scratch/transcripts.spades
    mkdir -p test/scratch/transcripts.spades
    cd test/scratch/transcripts.spades

    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config long_ssu_length=150 long_lsu_length=150 cleaning_protocol=assembly_no_ec --verbose > assembly.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -j 10 -s ../../../assembly.transcriptomic.snake -p --configfile ../../data/configs/spades.cdna.yaml --config long_ssu_length=150 long_lsu_length=150 cleaning_protocol=assembly_no_ec -n 2>&1 | grep '^Nothing to be done'"
    [ "$status" -eq 0 ]
    [ "${lines[0]:0:18}" == "Nothing to be done" ]
}
