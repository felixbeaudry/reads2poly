#!/bin/bash
#script to get dnds from outgroup
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/dnds_maker.sh loci.list allFasta nostop2buc.fasta
loc_list=$1
outDir=$2
outgroup_assembly=$3
outgroup_name=$4

if [#exists]; then
echo "Query exists"
else
##make database
/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/makeblastdb -in RB1_transcriptome_ref.fa -dbtype nucl
##make match file
/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/blastn -task megablast -query /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -db /ohta/felix.beaudry/assemblies/bucephalophorusTranscriptome/RB1_transcriptome_ref.fa -num_threads 10 -max_target_seqs 1 -outfmt 5 -dust no -gapopen 0  -gapextend 0 >nostop2buc/nostop2buc.out
fi

##parse
python blast_cleaner.py
#awk '$1 ~ "<Iteration_query-def>" ||  $1 ~ "<Hsp_hseq>" {print}' nostop2buc/nostop2buc.out | sed ':a;N;$!ba;s/<Hit_def>//g' | sed ':a;N;$!ba;s/<Iteration_query-def>//g' | sed ':a;N;$!ba;s/<Hsp_qseq>//g' | sed ':a;N;$!ba;s/<Hsp_hseq>//g' | sed ':a;N;$!ba;s/<\/Hsp_qseq>//g'| sed ':a;N;$!ba;s/<\/Iteration_query-def>//g'  | sed ':a;N;$!ba;s/      //g' | sed ':a;N;$!ba;s/<\/Hsp_hseq>\n//g' | sed ':a;N;$!ba;s/  Locus/\n\n>Locus/g' >nostop2buc.fasta


if [ -d "${outDir}/prank" ]; then
	rm ${outDir}/prank/*
else
	mkdir ${outDir}/prank
fi

samtools ${outgroup_assembly} ${loc}  | cat ${outDir}/${loc} > prank/${loc}

while read loc
do

##take output from fasta_maker and add outgroup
#samtools faidx nostop2roth.fasta ${loc} | cat $outdir/$loc >>${ind}/${ind}_${loc}.fasta

##add frame to outgroup sequence
perl /ohta/felix.beaudry/scripts/reads2poly/codoner.pl ${outDir}/${loc}.fasta >${outDir}/prank/${loc}_cod.fasta



done < $loc_list

##align to outgroup, in frame - using codon model -, with PRANK
/ohta/felix.beaudry/scripts/prank/bin/prank -d=${outDir}/prank/${loc}_cod.fasta -o=prank/${loc} -codon -F

##PAML

