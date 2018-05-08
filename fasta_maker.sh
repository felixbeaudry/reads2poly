#!/bin/bash
#script to make fastas of all sequences
#Felix Beaudry 8 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/fasta_maker.sh inds.list loci.list allFasta


ind_list=$1
loc_list=$2

outDir=$3

while read ind 
do
while read loc
do
rm $outDir/$loc.fasta
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i ${ind}/${loc}.fasta | cat >> $outDir/${loc}.fasta
#perl /ohta/felix.beaudry/scripts/reads2poly/codoner.pl
#align to outgroup, in frame, with PRANK
#send in to polymorphurama
done < $loc_list
done < $ind_list



