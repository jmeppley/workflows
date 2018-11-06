#!/usr/bin/env bats
setup() {
    which conda
    if [ "$?" -eq "1" ]; then
        unset conda
	export PATH=$_CONDA_ROOT:$PATH
    fi

    mkdir -p test/conda/envs
    ENV=nose
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ ! -e "$ENV_DIF" -o "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi

    source activate $ENV_DIR
}


@test "run nosetests" {
    run nosetests test/nose
    [ "$status" = 0 ]
}
