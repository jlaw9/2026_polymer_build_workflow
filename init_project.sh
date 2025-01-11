#!/bin/bash

# Parameters
runlevel=$1
n_sampled="${2:-6}"
datapath="${3:-"src/monomer_data/PolyID_master_unique.csv"}"
projname="${4:-"polyID"}"
outdir="${5:-"src"}"


# internal values
FLAG_TEST="--test"
FLAG_PROD="--production"

if [[ $runlevel == $FLAG_TEST ]]; then
    echo "Initializing test project..."
    projdir="${projname}_test"
    python -m src.init_signac -mdat $datapath -num $n_sampled --random --project-name $projdir -od $outdir
    cp -r 'src/templates' $projdir
    echo "Project directory '${projdir}' created for ${n_sampled} random chemistries"
elif [[ $runlevel == $FLAG_PROD ]]; then
    echo "Initializing production-scale project..."
    projdir="${projname}_production"
    python -m src.init_signac -mdat $datapath --project-name $projdir -od $outdir
    cp -r 'src/templates' $projdir
    echo "Project directory '${projdir}' created"
else
    echo "Invalid flag '$runlevel'; select either '$FLAG_TEST' or '$FLAG_PROD'"
fi
