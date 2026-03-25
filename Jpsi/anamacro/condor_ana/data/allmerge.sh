#!/usr/bin/env bash
set -euo pipefail


INDIR=result/data
OUTDIR=/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/data/all
OUTROOT=${OUTDIR}/mass_reco_run24pp_merged.root

hadd -f ${OUTROOT} ${INDIR}/*_mass_reco_merged.root
echo ${OUTROOT}
