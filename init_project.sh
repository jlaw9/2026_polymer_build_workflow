#!/bin/bash

# Parameters
runlevel=$1
n_sampled="${2:-6}"
datapath="${3:-"examples/monomers_example.csv"}"
projname="${4:-"polyID"}"
outdir="${5:-"."}"


# internal values
FLAG_TEST="--test"
FLAG_PROD="--production"

if [[ $runlevel == $FLAG_TEST ]]; then
    status="initialized"
    
    echo "Initializing test project..."
    projdir="${projname}_test"
    python -m src.init_signac -mdat $datapath -num $n_sampled --random --project-name $projdir -od $outdir || status="failed"
    
    if [ $status = "initialized" ]; then
        echo "Project directory '${projdir}' created for ${n_sampled} random chemistries"
    fi
elif [[ $runlevel == $FLAG_PROD ]]; then
    status="initialized"

    echo "Initializing production-scale project..."
    projdir="${projname}_production"
    python -m src.init_signac -mdat $datapath --project-name $projdir -od $outdir || status="failed"

    if [ $status = "initialized" ]; then
        echo "Project directory '${projdir}' created"
    fi
else
    echo "Invalid flag '$runlevel'; select either '$FLAG_TEST' or '$FLAG_PROD'"
fi
