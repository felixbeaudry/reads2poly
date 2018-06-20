#!/bin/bash
#script to make fastas of all sequences
#Felix Beaudry 8 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/fasta_maker.sh inds.list loci.list allFasta outgroup 

ind_list=$1
loc_list=$2
outDir=$3
outgroup=$4
#cutoff=$5

echo "Making $ind_list $loc_list fastas"

mkdir ${outDir}
awk 'split($1,a,"."){print a[1]}' ${loc_list} >${outDir}/${loc_list}
echo "Removing Previous Fasta Files"

while read loc
do
rm ${outDir}/${loc}.fasta
done < ${outDir}/${loc_list}

while read ind 
do
echo "Adding ${ind} sequences to fasta"
export PYTHONPATH=/usr/lib64/python2.7/site-packages/
while read loc
do
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i ${ind}/${loc}.fasta -c 60 2>${outDir}/errors.txt | cat >> ${outDir}/${loc}.fasta
#cat ${outDir}/${ind}/${loc}.fasta >> ${outDir}/${loc}.fasta
done < ${outDir}/${loc_list}
done < ${ind_list}

echo "Removing Outgroup Fasta Files"
while read outgr
do
while read loc
do
rm ${outDir}/${outgr}/${loc}.fasta
done < ${outDir}/${loc_list}
done < ${outgroup}

while read outgr
do
mkdir ${outDir}/${outgr}
while read loc
do
cat ${outgr}/${loc}.fasta ${outDir}/${loc}.fasta >${outDir}/${outgr}/${loc}.fasta
done < ${outDir}/${loc_list}
done < ${outgroup}

##take inds_list and print as one line with commas, and output into directory
#sed -E -e ':a;N;$!ba;s/\n/,/g' $ind_list >${outDir}/pop


