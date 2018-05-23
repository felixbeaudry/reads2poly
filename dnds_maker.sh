#!/bin/bash
#script to get dnds from outgroup
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/dnds_maker.sh auto_loci.list allFasta roth /ohta/felix.beaudry/assemblies/bucephalophorusTranscriptome/RB1_transcriptome_ref.fa
loc_list=$1
outDir=$2
outgroupName=$3
outgroupTranscriptome=$4

##make database
#/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/makeblastdb -in RB1_transcriptome_ref.fa -dbtype nucl
##make match file
mkdir ${outgroupName}
echo "Blasting ${outgroupName}"
/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/blastn -task megablast -query /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -db ${outgroupTranscriptome} -num_threads 10 -max_target_seqs 1 -outfmt 5 -dust no -gapopen 0 -gapextend 0 >${outgroupName}/${outgroupName}_blast.txt

##parse
echo "Parsing ${outgroupName}"
python blast_cleaner.py -i ${outgroupName}/${outgroupName}_blast.txt -o ${outgroupName}

mkdir ${outgroupName}/prank
while read loc
do
echo "Aligning ${loc}"
perl /ohta/felix.beaudry/scripts/reads2poly/codoner.pl ${outgroupName}/${loc}.fasta >${outgroupName}/${loc}_cod.fasta
samtools faidx ${outgroupName}/${loc}_cod.fasta ${outgroupName} | cat ${outDir}/${loc} > ${outgroupName}/${loc}.fasta
##align to outgroup, in frame - using codon model -, with PRANK
/ohta/felix.beaudry/scripts/prank/bin/prank -d=${outgroupName}/${loc}.fasta -o=${outgroupName}/prank/${loc} -codon -F
done < $loc_list


