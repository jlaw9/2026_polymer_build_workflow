#!/bin/bash

PROJDIR='polyID_test'

# Get current status of project
# python src/project.py -path $PROJDIR status --detailed --pretty

# Send all project jobs to the Cluster
python src/project.py -path $PROJDIR submit --parallel --partition blanca-shirts -o everything --job-output dumps.out # --pretend