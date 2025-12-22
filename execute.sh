#!/bin/bash

projdir="${1:-"polyID_test"}"
logfile="${2:-"polyid_build.out"}"

# Get current status of project
# python -m src.project -path $projdir status --detailed --pretty

# Send all project jobs to the Cluster
# python -m src.project -path $projdir submit --parallel --partition blanca-shirts -o everything --job-output dumps.out # --pretend # Shirts Blanca nodes
python -m src.project -path $projdir submit --parallel --bundle 6 -o everything --partition amilan --account ucb500_asc2 --job-output $logfile # Alpine allocation
