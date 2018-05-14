#!/bin/bash
#script to abba-baba files from single locus fastas
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/abba_maker.sh inds.list loci.list allFasta roth

ind_list=$1
loc_list=$2

outDir=$3
outgroup=$4

#check which loci have coverage in every individual

#take outgroup from dnds_maker

ind_count = (wc -l $ind_list)
if (wc -l ${loc}.fasta) >=  $ind_count 
	cat to abba.fasta


#cat all fasta that have complete coverage
