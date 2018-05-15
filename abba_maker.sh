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
rm ${ind}/${ind}_cat.fasta
done<$ind_list

ind_count = $(wc -l $loc_list | awk '{print $1*2 +1}')

while read loc
do

loc_ind_count = $(awk '$1 ~ ">" {print}' $outgroup/codon/prank/${loc}.fasta.best.fas | wc -l | awk '{print $1}')

##check which loci have coverage in every individual
if ["$loc_ind_count" -ge "$ind_count"]
then
	echo "Adding ${loc}"
	while read ind
	do
	samtools faidx roth/codon/prank/${loc}.fasta.best.fas ${ind} >>${ind}/${ind}_cat.fasta
	done <$ind_list
else 
	echo "${loc} not added"
fi

done <outgroup_loci.list


#while read ind
#do
# send to python script
	#add locus name
	#remove lines with >
	#remove end of line
##cat add individual to abba.fasta in order of sequences to $ind_list

#done <$ind_list

#run through the patD script (.r?)

