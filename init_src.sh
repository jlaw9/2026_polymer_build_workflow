#!/bin/sh

# initialize project-wide parameters and reaction templates
python -m src.parameters write -aow
python -m src.reactions
