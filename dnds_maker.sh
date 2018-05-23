#!/bin/bash
#script to get dnds from outgroup
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/dnds_maker.sh goodsex buc /ohta/felix.beaudry/assemblies/bucephalophorusTranscriptome/RB1_transcriptome_ref.fa

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/dnds_maker.sh goodsex roth /ohta/felix.beaudry/assemblies/rothTranscriptome/roth_female_ref_genewise.fa


outDir=$1
outgroupName=$2
outgroupTranscriptome=$3

##make database
#/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/makeblastdb -in RB1_transcriptome_ref.fa -dbtype nucl
##make match file
#mkdir ${outgroupName}
#echo "Blasting ${outgroupName}"
#/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/blastn -task megablast -query /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -db ${outgroupTranscriptome} -num_threads 10 -max_target_seqs 1 -outfmt 5 -dust no -gapopen 0 -gapextend 0 >${outgroupName}/${outgroupName}_blast.txt

##parse
echo "Parsing ${outgroupName}"
python /ohta/felix.beaudry/scripts/reads2poly/blast2fullfasta.py -i ${outgroupName}/${outgroupName}_blast.txt -o ${outgroupName}

ls ${outgroupName} >${outgroupName}/loc.list
#awk 'split($1,a,"."){print a[1]}' ${loc_list} >${outgroupName}/loc.list

while read loc
do
cat ${outDir}/${loc} >>${outgroupName}/${loc}
done < ${outgroupName}/loc.list

#while read loc
#do
#echo "Aligning ${loc}"
#perl /ohta/felix.beaudry/scripts/reads2poly/codoner.pl ${outgroupName}/${loc}.fasta >${outgroupName}/prank/${loc}.fasta
#perl /ohta/felix.beaudry/scripts/reads2poly/codoner.pl ${outDir}/${loc}.fasta  >> ${outgroupName}/prank/${loc}.fasta
##align to outgroup, in frame - using codon model -, with PRANK
#/ohta/felix.beaudry/scripts/prank/bin/prank -d=${outgroupName}/prank/${loc}.fasta -o=${outgroupName}/prank/${loc} -codon -F
#done < ${outgroupName}/loc.list


