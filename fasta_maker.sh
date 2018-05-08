#!/bin/bash
#script to make fastas of all sequences
#Felix Beaudry 8 May 2018

#run with: bash fasta_maker.sh inds.list loci.list outDir


ind_list=$1
loc_list=$2

outDir=$3
rm $outDir/*

while read ind 
do
while read loc
do
python /ohta/felix.beaudry/scripts/fasta_cleaner.py -i ${ind}/${loc}.fasta | cat >> $outDir/${loc}.fasta
done < $loc_list
done < $ind_list

#align to outgroup, in frame (?)

#send in to polymorphurama
