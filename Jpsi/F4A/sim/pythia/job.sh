#!/usr/bin/env bash

nevts="$1"
indst="$2"
outroot="$3"
startevt="$4"

echo "Nevts = ${nevts}"
echo "In-DST = ${indst}"
echo "Out-Root = ${outroot}"
echo "Start-event = ${startevt}"

#source /opt/sphenix/core/bin/sphenix_setup.csh -n new
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.516
# source /opt/sphenix/core/bin/sphenix_setup.sh -n new

# Additional commands for my local environment

# set workdir="/sphenix/user/hachiya/INTT/INTT/general_codes/hachiya/SiliconSeeding/SiliconSeedAna/"
# set workdir="/sphenix/user/mikemoto/SiCalo/siliconseedana"
myinstall="/sphenix/user/mikemoto/SiCalo/siliconseedana/install"
source /opt/sphenix/core/bin/setup_local.sh "${myinstall}"

echo "Fun4All_DST_SiliconSeedAna.C(${nevts},\"${indst}\",\"${outroot}\",${startevt})"
root -b -q "Fun4All_DST_SiliconSeedAna.C(${nevts},\"${indst}\",\"${outroot}\",${startevt})"

ls -lh "${outroot}" || echo "[WARN] output not found"

echo "all done"
date