#!/bin/sh

# initialize project-wide parameters and reaction templates
python src/parameters.py
python src/reactions.py

# create master monomer datafiles in standardized format
bash format_mdat.sh