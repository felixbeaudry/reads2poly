#!/bin/bash
#script to abba-baba files from single locus fastas
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/abba_maker.sh inds.list roth 

ind_list=$1
outgroup=$2

ls ${outgroup}/codon/prank/* | awk 'split($1,a,".") split(a[1],b,"/") {print b[4]}' >${outgroup}_loci.list

echo "Removing Past Files"
while read ind
do
rm ${ind}/${ind}_1_cat.fasta
rm ${ind}/${ind}_2_cat.fasta
done<$ind_list

rm ${outgroup}/${outgroup}_cat.fasta

while read loc
do

indcount=$(wc -l $ind_list | awk '{print $1*2 +1}')
locindcount=$(awk '$1 ~ ">" {print}' ${outgroup}/codon/prank/${loc}.fasta.best.fas | wc -l | awk '{print $1}')

##check which loci have coverage in every individual
if [[ $locindcount = $indcount ]]; then
echo "Adding ${loc}"
samtools faidx ${outgroup}/codon/prank/${loc}.fasta.best.fas ${outgroup} >>${outgroup}/${outgroup}_cat.fasta
while read ind
do
for hap in 1 2
do
samtools faidx ${outgroup}/codon/prank/${loc}.fasta.best.fas ${ind}_${hap} >>${ind}/${ind}_${hap}_cat.fasta
done
done <$ind_list
else 
echo "${loc} not added"
fi
done <${outgroup}_loci.list

rm sorted_abba.fasta
echo "Collapsing fastas"
while read ind
do
for hap in 1 2
do
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cat.py -i ${ind}/${ind}_${hap}_cat.fasta -n ${ind}_${hap} >>sorted_abba.fasta
done
done <$ind_list
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cat.py -i ${outgroup}/${outgroup}_cat.fasta -n ${outgroup} >>sorted_abba.fasta
#python /ohta/felix.beaudry/scripts/reads2poly/fasta_cat.py -i roth/roth_cat.fasta -n roth >>sorted_abba.fasta

samtools faidx sorted_abba.fasta NCROS7_2 FLJAS13_2 TXROS24_2 roth >abba_input.fasta

##rename sequences to P1, P2, P3, OUTGROUP
#sed ':a;N;$!ba;s/>NCROS7_2/>pop1/g' abba_input.fasta | sed ':a;N;$!ba;s/>FLJAS13_2/>pop2/g' 

#run through the patD script
Rscript --vanilla /ohta/felix.beaudry/scripts/reads2poly/patD.R abba_input.fasta >patD.txt

