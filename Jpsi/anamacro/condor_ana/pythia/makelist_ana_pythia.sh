#!/usr/bin/env bash

# RUN=${1:?usage: $0 RUNNUMBER}

INDIR=/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana/condor_merge
OUTDIR=/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana
OUTFILE=${OUTDIR}/input_pythia.list

: > "$OUTFILE" 
# mkdir -p ${OUTDIR}

ls -1 ${INDIR}/ana_merge_*.root | sort > ${OUTFILE}

echo ${OUTFILE}
