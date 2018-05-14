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
echo "Adding sequences"
perl /ohta/felix.beaudry/scripts/reads2poly/vcf2fasta_uni.pl  -v ${ind}.uni.vcf -o ${ind} -l ${ind}.logfile -a T 2>${ind}/errors.txt

while read loc
do
echo -e "${ind}\t${loc}"
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i ${ind}/${loc}.fasta 2>$outDir/errors.txt | cat >> $outDir/${loc}.fasta
#python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i 3_Index-10.TM4/Locus_9983_Transcript_2_Length_951_genewise_newlength_951_CDSlength_645.fasta | cat >> allFasta/Locus_9983_Transcript_2_Length_951_genewise_newlength_951_CDSlength_645.fasta 


done < $loc_list
done < $ind_list

#send in to polymorphurama
perl /ohta/felix.beaudry/scripts/reads2poly/polymorphurama_interpop.pl fasta $outDir


