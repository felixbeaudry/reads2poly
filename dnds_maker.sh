#!/bin/bash
#script to get dnds from outgroup
#Felix Beaudry 14 May 2018

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/dnds_maker.sh goodsex buc /ohta/felix.beaudry/assemblies/bucephalophorusTranscriptome/RB1_transcriptome_ref.fa 

#run with: bash /ohta/felix.beaudry/scripts/reads2poly/dnds_maker.sh goodsex roth /ohta/felix.beaudry/assemblies/rothTranscriptome/roth_female_ref_genewise.fa goodsex.fasta.list


outDir=$1
outgroupName=$2
outgroupTranscriptome=$3
locList=$4

mkdir ${outgroupName}
mkdir ${outgroupName}/perlocus
mkdir ${outgroupName}/overall

echo "Blasting ${loc}"
/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/blastn -task megablast -query /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -db ${outgroupTranscriptome} -num_threads 10 -max_target_seqs 1 -outfmt 5 -dust no -gapopen 0 -gapextend 0 >${outgroupName}/overall/${loc}.blast 2>${outgroupName}/overall/${loc}.blast.err
echo "Parsing ${loc}"
python /ohta/felix.beaudry/scripts/reads2poly/blast2fullfasta.py -i ${outgroupName}/overall/${loc}.blast -o ${outgroupName}/overall/${loc} 2>${outgroupName}/overall/${loc}.err

while read loc
do
##make database
#/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/makeblastdb -in ${outgroupTranscriptome} -dbtype nucl
echo "Blasting ${loc}"
/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/blastn -task megablast -query ${outDir}/${loc} -db ${outgroupTranscriptome} -num_threads 10 -max_target_seqs 1 -outfmt 5 -dust no -gapopen 0 -gapextend 0 >${outgroupName}/perlocus/${loc}.blast 2>${outgroupName}/perlocus/${loc}.blast.err
awk '$1 ~ "<Hit_def>" {print}' ${outgroupName}/perlocus/${loc}.blast | sort | uniq | awk 'split($1,a,">") split(a[2],b,"<") {print b[1]}' >${outgroupName}/perlocus/${loc}.blast.hits
done < $locList

#echo "Parsing ${loc}"
#python /ohta/felix.beaudry/scripts/reads2poly/blast2fullfasta.py -i ${outgroupName}/${loc}.blast -o ${outgroupName}/${loc}


#awk 'split($1,a,"."){print a[1]}' ${loc_list} >${outgroupName}/loc.list

#while read loc
#do
#cat ${outDir}/${loc} >>${outgroupName}/${loc}
#done < ${outgroupName}/loc.list




