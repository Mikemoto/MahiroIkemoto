#! /bin/bash

mode="jpsi"
# mode="e-"

output_ana="/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana_${mode}/merged_200k.root"
input_files_ana=""

# output_dst="/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/dst_jpsi/merged_10k.root"
# input_files_dst=""

while read line; do
  input_files_ana+="/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana_${mode}/condor/topo_ana_${line}.root "
  # input_files_dst+="/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/dst_jpsi/condor/DST_SiliconOnly_single_jpsi_vtxfixed_Caloevents_reconstructedInfo_${line}.root "
done < id_list.txt

hadd -f $output_ana $input_files_ana
# hadd -f $output_dst $input_files_dst

echo all done

date