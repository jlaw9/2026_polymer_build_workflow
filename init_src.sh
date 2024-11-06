#!/bin/sh

# initialize project-wide parameters and reaction templates
python -m src.parameters
python -m src.reactions

# create master monomer datafiles in standardized format
python -m src.format_data merge --glob "monomer_data_raw/*.csv" -rxns src/reactions/rxns_polyID.json -of PolyID_master_data.csv -od src/monomer_data -aow
python -m src.format_data merge -mdat "monomer_data_raw/nipu_urethanes.csv" -rxns src/reactions/rxns_polyID.json -of subsample_IPU.csv -od src/monomer_data -k 10 --random -aow