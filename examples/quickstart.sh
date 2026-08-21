#!/bin/bash
# End-to-end smoke test of the build workflow on 3 small chemistries.
#
#   bash examples/quickstart.sh          # build 3 chemistries as trimers
#   bash examples/quickstart.sh 5 10     # ...5 chemistries, DOP 10
#
# Expect a few minutes per chemistry: parameterization and packing dominate.
set -euo pipefail

N_CHEM="${1:-3}"
DOP="${2:-3}"
PROJ="quickstart_project"
MDAT="examples/monomers_example.csv"

echo "== 1/5 reaction templates =="
python -m src.reactions

echo "== 2/5 build parameters (DOP=$DOP) =="
python -m src.parameters write -aow --DOP "$DOP"

echo "== 3/5 format the monomer data =="
# adds the canonical/explicit SMILES columns that init_signac requires
python -m src.format_data sequence -mdat "$MDAT" -od examples -pf fmt -aow
MDAT_FMT="examples/monomers_example_fmt.csv"

echo "== 4/5 initialize project ($N_CHEM chemistries) =="
rm -rf "$PROJ"
python -m src.init_signac -mdat "$MDAT_FMT" -num "$N_CHEM" --project-name "$PROJ" -od .

echo "== 5/5 run the build =="
python -m src.project -path "$PROJ" -rxns src/reactions/rxns_polyID.json run -o everything

echo
echo "Done. Project status:"
python -m src.project -path "$PROJ" -rxns src/reactions/rxns_polyID.json status
echo
echo "MD inputs are under $PROJ/workspace/<job-id>/"
