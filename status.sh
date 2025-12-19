#!/bin/bash

projdir="${1:-"polyID_test"}"
rxnpath="${2:-"src/reactions/rxns_polyID.json"}"

python src/project.py \
    -path "$projdir" \
    -rxns "$rxnpath" \
    status \