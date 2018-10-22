
#perl /ohta/felix.beaudry/scripts/find_perfect.pl >find_perfect_list.vcf
#awk '{print $1}' find_perfect_list.vcf | sort | uniq >sex.list

while read ind
do
echo "sorting ${ind}"
samtools sort ${ind}_add.bam >${ind}_sort.bam
java -jar /ohta/apps/GenomeAnalysisTK.jar -T SelectVariants -R /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -V ${ind}.uni.vcf -o ${ind}.sex.vcf -L sex.list
echo "phasing ${ind}"
/ohta/apps/hapcut/extractHAIRS --VCF ${ind}.sex.vcf --bam ${ind}_sort.bam --maxIS 400 --ref /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa > ${ind}.fragmatrix
/ohta/apps/hapcut/HAPCUT --fragments ${ind}.fragmatrix --VCF ${ind}.sex.vcf --output ${ind}.hapcut --maxmem 30000 > ${ind}.uni_hapcut.log
perl /ohta/felix.beaudry/scripts/extract_seq_perfect_sites.pl find_perfect_list.vcf ${ind}.hapcut ${ind}.sex.vcf >${ind}.fasta
done <rna.inds
