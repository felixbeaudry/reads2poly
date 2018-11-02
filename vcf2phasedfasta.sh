

perl /ohta/felix.beaudry/scripts/find_perfect.pl >find_perfect_list.vcf
awk '{print $1}' find_perfect_list.vcf | sort | uniq >sex.list

rm sex.fasta.err
while read ind
do
echo "Variant calling for ${ind}"
java -Djava.io.tmpdir=tmp -jar /ohta/apps/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -I ${ind}_ddpl.bam -o ${ind}.uni.vcf --output_mode EMIT_ALL_SITES -nct 10 --genotyping_mode DISCOVERY

#samtools sort ${ind}_add.bam >${ind}_sort.bam

java -jar /ohta/apps/GenomeAnalysisTK.jar -T SelectVariants -R /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -V ${ind}.uni.vcf -o ${ind}.sex.vcf -L sex.list --selectTypeToExclude INDEL

echo "Phasing ${ind}"
/ohta/apps/hapcut/extractHAIRS --VCF ${ind}.sex.vcf --bam ${ind}_sort.bam --maxIS 400 --ref /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa > ${ind}.fragmatrix
/ohta/apps/hapcut/HAPCUT --fragments ${ind}.fragmatrix --VCF ${ind}.sex.vcf --output ${ind}.hapcut --maxmem 30000 > ${ind}.uni_hapcut.log

perl /ohta/felix.beaudry/scripts/extract_seq_perfect_sites.pl find_perfect_list.vcf ${ind}.hapcut ${ind}.sex.vcf >${ind}.sex.fasta 2>${ind}.sex.fasta.err

perl /ohta/felix.beaudry/scripts/reads2poly/fasterLiner.pl -i ${ind}.sex.fasta

#mkdir phase/${ind}
rm phase/${ind}/*; 
for chrom in X Y
do

echo "Fasting ${ind} ${chrom}"	
while read loc 
do
samtools faidx ${ind}.sex.fasta.sort ${loc}_${chrom} 2>>phase/${ind}/${ind}.fai.out | cat >phase/${ind}/${loc}_${chrom}.fasta
cat phase/${ind}/${loc}_${chrom}.fasta | grep -v '^>' | grep '^.' | tr -d '[:blank:]' | tr n - | tr , - > phase/${ind}/${loc}_${chrom}.noname.fasta
echo ">${ind}_${chrom}chrom" | cat - phase/${ind}/${loc}_${chrom}.noname.fasta > temp && mv temp phase/${ind}/${loc}_${chrom}.fasta
rm phase/${ind}/${loc}_${chrom}.noname.fasta
done < sex.list
done
done < male.inds


for chrom in X Y
do
while read loc
do 
for outgroup in bucephalophorus rothschildianus
rm subsets/rna/${chrom}phase/${loc}.fasta
cat /ohta/felix.beaudry/fastas/bucephalophorus/${loc}.fasta >> subsets/rna/${chrom}phase/${loc}.fasta
cat /ohta/felix.beaudry/fastas/rothschildianus/${loc}.fasta >> subsets/rna/${chrom}phase/${loc}.fasta
while read ind
do
python /ohta/felix.beaudry/scripts/reads2poly/fasta_cleaner.py -i phase/${ind}/${loc}_${chrom}.fasta -c 60 >> subsets/rna/${chrom}phase/${loc}.fasta
done <male.inds
done <sex.list
done


for loc in X Y
do
for pop in pop
do
for subset in rna
do
mkdir subsets/${subset}/${loc}phase
cp ${pop}.pop subsets/${subset}/${loc}phase/${pop}.pop
for outgroup in rothschildianus bucephalophorus
do
perl /ohta/felix.beaudry/scripts/reads2poly/polymorphurama_interpop.pl -i subsets/${subset}/${loc}phase/ -p ${pop}.pop -S ${subset} -o ${outgroup} -C ${loc} -c ${loc}
done
done
done
done


