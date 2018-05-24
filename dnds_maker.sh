#!/bin/bash
#script to get dnds from outgroup
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/dnds_maker.sh goodsex buc /ohta/felix.beaudry/assemblies/bucephalophorusTranscriptome/RB1_transcriptome_ref.fa goodsex.fasta.list

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/dnds_maker.sh goodsex roth /ohta/felix.beaudry/assemblies/rothTranscriptome/roth_female_ref_genewise.fa goodsex.fasta.list

outDir=$1
outgroupName=$2
outgroupTranscriptome=$3
locList=$4

##make database
#/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/makeblastdb -in ${outgroupTranscriptome} -dbtype nucl

mkdir ${outgroupName}
echo "Blasting ${outgroupName}"
/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/blastn -task megablast -query /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -db ${outgroupTranscriptome} -num_threads 10 -max_target_seqs 1 -outfmt 5 -dust no -gapopen 0 -gapextend 0 >${outgroupName}/${outgroupName}.blast 2>${outgroupName}/${outgroupName}.blast.err
echo "Parsing ${outgroupName}"
python /ohta/felix.beaudry/scripts/reads2poly/blast2fullfasta.py -i ${outgroupName}/${outgroupName}.blast -o ${outgroupName} -s ${outgroupName} 2>${outgroupName}/${outgroupName}.err

awk 'split($1,a,"."){print a[1]}' ${locList} >${outgroupName}/loc.list

mkdir ${outDir}/${outgroupName}
while read loc
do
cat ${outgroupName}/${loc}.fasta ${outgroupName}/${loc}.fasta >>${outDir}/${outgroupName}/${loc}.fasta
done < ${outgroupName}/loc.list






