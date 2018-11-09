setup() {
    which conda
    if [ "$?" -eq "1" ]; then
        unset conda
	export PATH=$_CONDA_ROOT:$PATH
    fi

    mkdir -p test/conda/envs
    ENV=illumina.qc
    ENV_DIR=`pwd`/test/conda/envs/$ENV
    ENV_FILE=test/conda/${ENV}.yml
    if [ "$ENV_FILE" -nt "$ENV_DIR" ]; then
        rm -rf $ENV_DIR
        conda env create -f $ENV_FILE -p $ENV_DIR --force --quiet > test/conda/envs/.create.$ENV 2>&1
    fi
    
    source activate $ENV_DIR
}

@test "Find top hits in a database" {
    rm -rf test/scratch/top
    mkdir -p test/scratch/top
    cd test/scratch/top
    run bash -c "snakemake -s ../../../annotation.tophits.snake --configfile ../../data/configs/reads.top.RS.yaml -p -k -j 5 > top.annot.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.tophits.snake --configfile ../../data/configs/reads.top.RS.yaml -p -k -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Check tophit outputs" {
    run bash -c "diff <(sort test/scratch/top/ALOHA_XVII_1_04.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts) test/bats/outputs/top/ALOHA_XVII_1_04.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts.sort"
    [ "$status" -eq 0 ]
    run bash -c "diff <(sort test/scratch/top/ALOHA_XVII_1_15.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts) test/bats/outputs/top/ALOHA_XVII_1_15.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts.sort"
    [ "$status" -eq 0 ]
    run bash -c "diff <(sort test/scratch/top/ALOHA_XVII_1_16.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts) test/bats/outputs/top/ALOHA_XVII_1_16.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts.sort"
    [ "$status" -eq 0 ]
    run bash -c "diff <(sort test/scratch/top/ALOHA_XVII_1_18.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts) test/bats/outputs/top/ALOHA_XVII_1_18.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts.sort"
    [ "$status" -eq 0 ]
    run bash -c "diff <(sort test/scratch/top/ALOHA_XVII_1_21.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts) test/bats/outputs/top/ALOHA_XVII_1_21.vs.RefSeq.lastx._E0.01_F0_I90.tophit.all.hitid.counts.sort"
    [ "$status" -eq 0 ]
}

@test "Find top hits in a database with QC" {
    rm -rf test/scratch/topqc
    mkdir -p test/scratch/topqc
    cd test/scratch/topqc
    run bash -c "snakemake -s ../../../annotation.tophits.snake --configfile ../../data/configs/reads.top.QC.RS.yaml -p -k -j 5 > top.annot.log 2>&1"
    [ "$status" -eq 0 ]
    run bash -c "snakemake -s ../../../annotation.tophits.snake --configfile ../../data/configs/reads.top.QC.RS.yaml -p -k -n 2>&1 | grep '^Nothing'"
    [ "$status" -eq 0 ]
    [ "${lines[0]}" == "Nothing to be done." ]
}

@test "Check tophit QC outputs" {
    run bash -c "diff <(sort test/scratch/topqc/ALOHA_XVII_1_04_DNA.vs.RefSeq.lastx._B50_F0.tophit.all.hitid.counts) test/bats/outputs/topqc/ALOHA_XVII_1_04_DNA.vs.RefSeq.lastx._B50_F0.tophit.all.hitid.counts.sort | grep -c '^>'"
    [ "${lines[0]}" -lt 75 ]
    run bash -c "diff <(sort test/scratch/topqc/ALOHA_XVII_1_04_DNA.vs.RefSeq.lastx._B50_F0.tophit.all.hitid.counts) test/bats/outputs/topqc/ALOHA_XVII_1_04_DNA.vs.RefSeq.lastx._B50_F0.tophit.all.hitid.counts.sort | grep -c '^<'"
    [ "${lines[0]}" -lt 75 ]
    run bash -c "diff <(sort test/scratch/topqc/ALOHA_XVII_1_15_DNA.vs.RefSeq.lastx._B50_F0.tophit.all.hitid.counts) test/bats/outputs/topqc/ALOHA_XVII_1_15_DNA.vs.RefSeq.lastx._B50_F0.tophit.all.hitid.counts.sort | grep -c '^>'"
    [ "${lines[0]}" -lt 75 ]
    run bash -c "diff <(sort test/scratch/topqc/ALOHA_XVII_1_15_DNA.vs.RefSeq.lastx._B50_F0.tophit.all.hitid.counts) test/bats/outputs/topqc/ALOHA_XVII_1_15_DNA.vs.RefSeq.lastx._B50_F0.tophit.all.hitid.counts.sort | grep -c '^<'"
    [ "${lines[0]}" -lt 75 ]
}


