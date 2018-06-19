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

export PYTHONPATH=/usr/lib64/python2.7/site-packages/
while read ind 
do
echo "Extracting ${ind} sequences from vcf"
perl /ohta/felix.beaudry/scripts/reads2poly/vcf2fasta_uni.pl  -v /ohta/felix.beaudry/alignments/NCF1_nostop/RNA/${ind}.uni.vcf -o ${outDir}/${ind} -l ${outDir}/${ind}.logfile -a T -m 60 -d 10 -g 30 2>${outDir}/${ind}_vcferrors.txt
echo "Adding ${ind} sequences to fasta"
while read loc
do
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i ${outDir}/${ind}/${loc}.fasta -c 60 2>${outDir}/errors.txt | cat >> ${outDir}/${loc}.fasta
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


