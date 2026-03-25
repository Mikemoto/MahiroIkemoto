#!/usr/bin/env bash

in_root="$1"

echo "input file = ${in_root}"

# source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.516
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.527
# source /opt/sphenix/core/bin/sphenix_setup.sh -n new

# Additional commands for my local environment

echo "DrawMassDis_data.C(\"${in_root}\")"
root -b -q "DrawMassDis_data.C(\"${in_root}\",1)"

echo "all done"
date
