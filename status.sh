#!/bin/bash

projdir="${1:-"quickstart_project"}"
rxnpath="${2:-"src/reactions/rxns_polyID.json"}"

python -m src.project \
    -path $projdir \
    -rxns $rxnpath \
    status \