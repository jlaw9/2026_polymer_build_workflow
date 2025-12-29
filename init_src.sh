#!/bin/sh

# initialize project-wide parameters and reaction templates
python -m src.parameters write -aow
python -m src.reactions

# create master monomer datafiles in standardized format
bash format_mdat.sh