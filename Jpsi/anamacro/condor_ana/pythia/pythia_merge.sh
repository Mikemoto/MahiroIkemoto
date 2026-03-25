
#!/usr/bin/env bash
set -euo pipefail

INDIR=result/sim/ana_condor
OUTDIR=/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/sim
OUTROOT=${OUTDIR}/pythia_mass_reco_100m_merged.root

hadd -f ${OUTROOT} ${INDIR}/ana_merge_*.root
echo ${OUTROOT}
