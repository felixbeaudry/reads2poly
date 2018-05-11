#!/bin/bash
#script to make fastas of all sequences
#Felix Beaudry 8 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/fasta_maker.sh inds.list loci.list allFasta


ind_list=$1
loc_list=$2

outDir=$3

echo "Removing Previous Fasta Files"
while read loc
do
rm $outDir/$loc.fasta
done < $loc_list

while read ind 
do
echo "Adding sequences from $ind"
while read loc
do
#if seq exists
	python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i ${ind}/${loc}.fasta 2>$outDir/errors.txt | cat >> $outDir/${loc}.fasta
# if not, skip


#find outgroup sequence using BLAST
#add frame to outgroup sequence
#perl /ohta/felix.beaudry/scripts/reads2poly/codoner.pl $outDir/${loc}.fasta
#align to outgroup, in frame - using codon model -, with PRANK
#/ohta/felix.beaudry/scripts/prank/bin/prank -d=.fasta -o=prank/.fasta -codon -F
#send in to polymorphurama
#perl /ohta/felix.beaudry/scripts/polymorphurama_FM_xyy.pl fasta allFasta/
#send into PAML
done < $loc_list
done < $ind_list



