#!/bin/bash

projdir="${1:-"polyID_test"}"
rxnpath="${2:-"src/reactions/rxns_polyID.json"}"

python -m src.project \
    -path "$projdir" \
    -rxns "$rxnpath" \
    status \