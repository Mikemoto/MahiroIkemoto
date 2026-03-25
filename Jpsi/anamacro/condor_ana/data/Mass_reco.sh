#!/usr/bin/env bash
set -euo pipefail

RUN=${1:?usage: $0 RUNNUMBER}

MACRODIR=/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro
JOBFILE=${MACRODIR}/mass_ana_condor.job

OUTDIR=${MACRODIR}/result/ana_condor
LOGDIR=${OUTDIR}/log
mkdir -p ${LOGDIR}

WAITLOG=${LOGDIR}/wait_${RUN}_$$_$(date +%s).log

# (1) makelist
${MACRODIR}/makelist_ana.sh ${RUN}

# (2) submit: 
TMPJOB=$(mktemp /tmp/mass_ana_${RUN}_XXXX.job)
sed "s@^Log *=.*@Log = ${WAITLOG}@g" ${JOBFILE} > ${TMPJOB}
SUBMIT_OUT=$(condor_submit -append "RUN=${RUN}" ${TMPJOB})
echo "${SUBMIT_OUT}"

CLUSTER=$(echo "${SUBMIT_OUT}" | awk '/submitted to cluster/ {print $NF}' | tr -d '.')
echo "cluster=${CLUSTER}"
condor_q ${CLUSTER}

# (3) wait
sleep 2
condor_wait ${WAITLOG}

# (4) merge
${MACRODIR}/merge.sh ${RUN}

# (5) draw
cd ${MACRODIR}
root -q -b "DrawHist.C(${RUN})"
