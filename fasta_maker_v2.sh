#!/bin/bash
#script to make fastas of all sequences
#Felix Beaudry 8 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/fasta_maker.sh inds.list loci.list allFasta outgroup 80

ind_list=$1
loc_list=$2
outDir=$3
outgroup=$4
cutoff=$5

echo "Making $ind_list $loc_list fastas"

mkdir ${outDir}
awk 'split($1,a,"."){print a[1]}' ${loc_list} >${outDir}/${loc_list}
echo "Removing Previous Fasta Files"

while read loc
do
rm ${outDir}/${loc}.fasta
done < ${outDir}/${loc_list}

cat ${outgroup} ${ind_list} > ${outDir}/${ind_list}

while read ind 
do
echo "Adding ${ind} sequences to fasta"
export PYTHONPATH=/usr/lib64/python2.7/site-packages/
while read loc
do
#python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i ${ind}/${loc}.fasta -c ${cutoff} 2>${outDir}/errors.txt | cat >> ${outDir}/${loc}.fasta
cat ${ind}/${loc}.fasta >> ${outDir}/${loc}.fasta
done < ${outDir}/${loc_list}
done < ${outDir}/${ind_list}



