#!/bin/bash
#script to abba-baba files from single locus fastas
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/abba_maker.sh inds.list roth abba_file.txt

ind_list=$1
outgroup=$2
abba_file=$3

ls ${outgroup}/codon/prank/*.fas | awk 'split($1,a,".") split(a[1],b,"/") {print b[4]}' >${outgroup}_loci.list
#ls roth/codon/prank/*.fas | awk 'split($1,a,".") split(a[1],b,"/") {print b[4]}' >roth_loci.list


echo "Removing Past Files"
rm ${outgroup}/${outgroup}_cat.fasta
rm ${outgroup}/${outgroup}_collapse.fasta
while read ind
do
for hap in 1 2
do
rm ${ind}/${ind}_${hap}_cat.fasta
rm ${ind}/${ind}_${hap}_collapse.fasta
done
done<$ind_list

##check which loci have coverage in every individual
echo "Adding Loci"
while read loc
do
indcount=$(wc -l $ind_list | awk '{print $1*2 +1}')
locindcount=$(awk '$1 ~ ">" {print}' ${outgroup}/codon/prank/${loc}.fasta.best.fas | wc -l | awk '{print $1}')
if [[ $locindcount = $indcount ]]; then
echo "Adding ${loc}"
samtools faidx ${outgroup}/codon/prank/${loc}.fasta.best.fas ${outgroup} >>${outgroup}/${outgroup}_cat.fasta
while read ind
do
for hap in 1 2
do
samtools faidx ${outgroup}/codon/prank/${loc}.fasta.best.fas ${ind}_${hap} >>${ind}/${ind}_${hap}_cat.fasta
done
done <${ind_list}
#else 
#echo "${loc} not added"
fi
done <${outgroup}_loci.list

echo "Collapsing fastas"
while read ind
do
for hap in 1 2
do
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cat.py -i ${ind}/${ind}_${hap}_cat.fasta -n ${ind}_${hap} >${ind}/${ind}_${hap}_collapse.fasta
done
done <$ind_list
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cat.py -i ${outgroup}/${outgroup}_cat.fasta -n ${outgroup} >${outgroup}/${outgroup}_collapse.fasta

rm abba_input.fasta
##(((P1, P2), P3), OUTGROUP
while read abba
do
samtools faidx ${abba}/${abba}_2_collapse.fasta ${abba}_2 >>abba_input.fasta
done <$abba_file
samtools faidx ${outgroup}/${outgroup}_collapse.fasta ${outgroup} >>abba_input.fasta

#run through the patD script
echo "Calculating Patterson's D"
Rscript --vanilla /ohta/felix.beaudry/scripts/reads2poly/patd.R abba_input.fasta >patD.txt

