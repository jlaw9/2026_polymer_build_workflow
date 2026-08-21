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

echo "== 1/4 reaction templates =="
python -m src.reactions

echo "== 2/4 build parameters (DOP=$DOP) =="
python -m src.parameters write -aow --DOP "$DOP"

echo "== 3/4 initialize project ($N_CHEM chemistries) =="
rm -rf "$PROJ"
python -m src.init_signac -mdat "$MDAT" -num "$N_CHEM" --project-name "$PROJ" -od .

echo "== 4/4 run the build =="
python -m src.project -path "$PROJ" -rxns src/reactions/rxns_polyID.json run -o everything

echo
echo "Done. Project status:"
python -m src.project -path "$PROJ" -rxns src/reactions/rxns_polyID.json status
echo
echo "MD inputs are under $PROJ/workspace/<job-id>/"
