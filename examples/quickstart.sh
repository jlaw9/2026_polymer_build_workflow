#!/bin/bash
# End-to-end run of the build workflow on example chemistries.
#
#   bash examples/quickstart.sh               # build 3 chemistries as trimers
#   bash examples/quickstart.sh 5 10          # ...5 chemistries, DOP 10
#   bash examples/quickstart.sh 5 10 10000    # ...10k atoms per system
#
# Expect a few minutes per chemistry: parameterization and packing dominate.
set -euo pipefail

N_CHEM="${1:-3}"
DOP="${2:-3}"
N_ATOMS_MAX="${3:-5000}"
PROJ="quickstart_project"
MDAT="examples/monomers_example.csv"
RXNS="src/reactions/rxns_polyID.json"
PARAMS="parameters_quickstart"

echo "== 1/5 initialize reaction templates =="
python -m src.reactions

echo "== 2/5 build parameters (DOP=$DOP) =="
python -m src.parameters write -aow \
    --DOP "$DOP" \
    --n-atoms-max "$N_ATOMS_MAX" \
    -namdat "$PARAMS"

echo "== 3/5 format the monomer data =="
# inject canonical and explicit SMILES columns required to specify explicit chemistry 
python -m src.format_data sequence -aow \
    -mdat "$MDAT" \
    -od examples \
    -pf fmt
MDAT_FMT="examples/monomers_example_fmt.csv"

echo "== 4/5 initialize project '$PROJ' with $N_CHEM chemistries =="
rm -rf "$PROJ" # remove any prior project with the same name
python -m src.init_signac -od . \
    --project-name "$PROJ" \
    -mdat "$MDAT_FMT" \
    -num "$N_CHEM" \
    --parameters-swept-name "$PARAMS"

echo "== 5/5 run local system build =="
python -m src.project \
    -path "$PROJ" \
    -rxns "$RXNS" \
    run -o everything

echo
echo "Done. Project status:"
python -m src.project -path "$PROJ" -rxns "$RXNS" status
echo
echo "MD inputs are under $PROJ/workspace/<job-id>/"
