#!/usr/bin/env bash
set -euo pipefail

RUN=${1:?usage: $0 RUNNUMBER}

INDIR=result/ana_condor/${RUN}
OUTDIR=/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/data
OUTROOT=${OUTDIR}/${RUN}_mass_reco_merged.root

hadd -f ${OUTROOT} ${INDIR}/ana_${RUN}_*.root
echo ${OUTROOT}
