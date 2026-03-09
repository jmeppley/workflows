#!/usr/bin/env bats
setup() {
    eval "$(conda shell.bash hook)"
    mkdir -p test/conda/envs
    ENV=nose
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ ! -e "$ENV_DIF" -o "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --quiet > test/conda/envs/.create.$ENV 2>&1
    fi

    conda activate $ENV_DIR
}


@test "run nosetests" {
    mkdir -p test/scratch
    run nosetests test/nose > test/scratch/nose.log 2>&1
    [ "$status" = 0 ]
}
