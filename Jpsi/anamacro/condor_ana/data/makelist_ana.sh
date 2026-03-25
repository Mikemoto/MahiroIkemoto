#!/usr/bin/env bash

RUN=${1:?usage: $0 RUNNUMBER}

INDIR=/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/makedata/condor/data/${RUN}
OUTDIR=/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/newgeo/calocalib/recalc_charge/list
OUTFILE=${OUTDIR}/input_${RUN}.list

mkdir -p ${OUTDIR}

ls -1 ${INDIR}/ana_${RUN}_*.root | sort > ${OUTFILE}

echo ${OUTFILE}
