setup() {
    mkdir -p test/conda/envs
    ENV=annotate
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi
    source activate $ENV_DIR
    conda install -y -c bioconda -c conda-forge pytest
}

@test "Pytests in python code" {
    python -m pytest --pyargs python
    [ "$status" -eq 0 ]
}

