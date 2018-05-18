#!/bin/bash
#script to make fastas of all sequences
#Felix Beaudry 8 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/fasta_maker.sh inds.list loci.list allFasta cutoff

ind_list=$1
loc_list=$2
outDir=$3
cutoff=$4

mkdir ${outDir}
awk 'split($1,a,"."){print a[1]}' ${loc_list} >${outDir}/${loc_list}
echo "Removing Previous Fasta Files"
while read loc
do
rm $outDir/$loc.fasta
done < ${outDir}/${loc_list}

while read ind 
do
echo "Adding sequences"
if [ -d "${ind}" ]; then
echo -e "${ind}\tfastas exist"
else
perl /ohta/felix.beaudry/scripts/reads2poly/vcf2fasta_uni.pl  -v ${ind}.uni.vcf -o ${ind} -l ${ind}.logfile -a T 2>${ind}_vcferrors.txt
fi
while read loc
do
echo -e "${ind}\t${loc}"
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i ${ind}/${loc}.fasta -c ${cutoff} 2>$outDir/errors.txt | cat >> $outDir/${loc}.fasta

done < ${outDir}/${loc_list}
done < $ind_list

##take inds_list and print as one line with commas, and output into directory
sed -E -e ':a;N;$!ba;s/\n/,/g' $ind_list >${outDir}/pop
##send in to polymorphurama
perl /ohta/felix.beaudry/scripts/reads2poly/polymorphurama_interpop.pl fasta ${outDir} pop 


