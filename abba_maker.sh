#!/bin/bash
#script to abba-baba files from single locus fastas
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/abba_maker.sh inds.list roth 

ind_list=$1

outgroup=$2

ls $outgroup/codon/prank/* | awk 'split($1,a,".") split(a[1],b,"/") {print b[4]}' >outgroup_loci.list

echo "Removing Past Files"
while read ind
do
rm ${ind}/${ind}_1_cat.fasta
rm ${ind}/${ind}_2_cat.fasta
done<$ind_list

while read loc
do

indcount=$(wc -l $ind_list | awk '{print $1*2 +1}')
locindcount=$(awk '$1 ~ ">" {print}' $outgroup/codon/prank/${loc}.fasta.best.fas | wc -l | awk '{print $1}')

##check which loci have coverage in every individual
if [[ $locindcount = $indcount ]]; then
echo "Adding ${loc}"
while read ind
do
for hap in 1 2
do
samtools faidx roth/codon/prank/${loc}.fasta.best.fas ${ind}_${hap} >>${ind}/${ind}_${hap}_cat.fasta
done
done <$ind_list
else 
echo "${loc} not added"
fi

done <outgroup_loci.list

rm abba.fasta
while read ind
do
for hap in 1 2
do
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cat.py -i ${ind}/${ind}_${hap}_cat.fasta -n ${ind}_${hap} >>sorted_abba.fasta
done
done <$ind_list

#rename sequences to outgroup pop1 pop2 pop3

#Rscript --vanilla patD.R abba.fasta >patD.txt
#run through the patD script (.r?)

