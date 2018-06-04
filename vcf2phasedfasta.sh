
perl /ohta/felix.beaudry/scripts/find_perfect.pl >find_perfect_list.vcf


while read ind
do

#filt
#clean
samtools sort ${ind}_add.bam >${ind}_addsort.bam
/ohta/apps/hapcut/extractHAIRS --VCF ${ind}.clean.vcf --bam ${ind}_addsort.bam --maxIS 400 --ref /ohta/felix.beaudry/transcriptomeassembly/NCF1_a.fa --mbq 20 --mmq 40 > ${ind}.fragmatrix
/ohta/apps/hapcut/HAPCUT --fragments ${ind}.full.fragmatrix --VCF ${ind}.clean.vcf --output ${ind}.full.hapcut --maxiter 100 --maxmem 30000 --mbq 20 --qvoffset 60 > ${ind}.uni_hapcut.log
done
perl /ohta/felix.beaudry/scripts/extract_seq_perfect_sites.pl /ohta/felix.beaudry/alignments/NCF1_a/perfect_sites_aa.vcf ${ind}.uni.hapcut ${ind}.uni.vcf >${ind}.fasta

perl extract_seq_perfect_sites_swap.pl /ohta/felix.beaudry/beyond/popProject/find_perfect_list.vcf /ohta/felix.beaudry/alignments/NCF1/pop/${ind}.uni.hapcut /ohta/felix.beaudry/alignments/NCF1/pop/${ind}.uni.vcf >/ohta/felix.beaudry/beyond/popProject/idealsites/${ind}_pollen.fasta 2>${ind}.err

