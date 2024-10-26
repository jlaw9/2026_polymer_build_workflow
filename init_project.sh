#!/bin/bash

# Parameters
runlevel=$1
projname="${2:-polyID}"

DATAPATH="src/monomer_data/PolyID_master_data.csv" # hard-coded for now
OUTDIR="src" # hard-coded for now


# internal values
FLAG_TEST="--test"
FLAG_PROD="--production"

if [[ $runlevel == $FLAG_TEST ]]; then
    echo "Initializing test project..."
    projdir="${projname}_test"
    python -m src.init_signac -mdat $DATAPATH -num 12 --random --project-name $projdir -od $OUTDIR
    cp -r 'src/templates' $projdir
    echo "Project directory '${projdir}' created"
elif [[ $runlevel == $FLAG_PROD ]]; then
    echo "Initializing production-scale project..."
    projdir="${projname}_production"
    python -m src.init_signac -mdat $DATAPATH --project-name $projdir -od $OUTDIR
    cp -r 'src/templates' $projdir
    echo "Project directory '${projdir}' created"
else
    echo "Invalid flag '$runlevel'; select either '$FLAG_TEST' or '$FLAG_PROD'"
fi
