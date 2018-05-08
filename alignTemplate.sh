
#a quick little script to get read of really bad reads
perl /ohta/felix.beaudry/scripts/clean_B_N_HiSeq.pl -p T -f "*_R1.fastq"
#stats for your reads, use this step to decide clip length for TrimGalore
/ohta/felix.beaudry/scripts/FastQC/fastqc ${ind}_clean.fastq.gz -t 10
for ind in
do
#Trim out adaptors ++
/ohta/felix.beaudry/scripts/FastQC/TrimGalore-0.4.3/trim_galore --paired -o ./trimmed/ --retain_unpaired --clip_r1 10 --three_prime_clip_r1 10  --clip_r2 10 --three_prime_clip_r2 10  ${ind}_R1.fastq.gz ${ind}_R2.fastq.gz


# this is the best practices recommendation for genomic alignment. Each of the two specialized in a different way to align sequences (one is fast, the second is accurate), running them in sequence saves time but also allows for accuracy
/ohta/apps/bwa-0.7.15/bwa mem -t 10 -M assembly.fasta _R1.fastq _R2.fastq | samtools view -Sb -o _bwa.bam  ##-M Mark shorter split hits as secondary
#useful for divergent reads 
python /ohta/apps/stampy-1.0.30/stampy.py -g -h -M _bwa.bam -o _stampy.sam --bamkeepgoodreads --bwamark -t 10

#this is the best practices recommendation for transcript/RNAseq alignment. It specialized in considering alternative splicing and the removal of introns.
/ohta/apps/STAR-master/bin/Linux_x86_64/STAR --genomeDir --readFilesIn .fastq --outFileNamePrefix --runThreadN 10

#add information about the read group (eg what, where, when, by whom)
java -jar /ohta/apps/picard/build/libs/picard.jar AddOrReplaceReadGroups I=.sam O=.bam  RGLB=  RGPU= RGSM= RGPL=Illumina SO=coordinate # TMP_DIR= # LB library (eg project name) PU platform (eg sequencing location) SM sample (eg individual)
#mark reads that look suspiciously similar
java -jar /ohta/apps/picard/build/libs/picard.jar MarkDuplicates I=.bam I=.bam O=.bam M=.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT  VALIDATION_STRINGENCY=SILENT TMP_DIR=/ohta/felix.beaudry/tmp/ MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

# these steps give you the number of reads aligned for each locus; the first makes an index of file, the second summarizes it in a txt file
samtools index .bam 
samtools idxstats .bam >.txt

# cleaning step, specific to RNAseq, again specializing in removing alternative splice junctions
java -jar /ohta/apps/GenomeAnalysisTK.jar -T SplitNCigarReads -R .fasta -I .bam -o .bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 


# for individuals, use this. this makes a file calling all the SNPs relative to the reference 
java -Djava.io.tmpdir=tmp -jar /ohta/apps/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa  -I ${ind}_bwa_add.bam -o ${ind}.uni.vcf  --output_mode EMIT_ALL_SITES

done

# for many individuals, use this. this makes a file calling all the SNPs relative to the reference 
java -Djava.io.tmpdir=tmp -jar /ohta/apps/GenomeAnalysisTK.jar -T HaplotypeCaller -R .fasta -I .bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o .g.vcf -ERC GVCF 
#this makes a file containing a concatination of the .vcf s made in the previous step, so that you can compare between samples. 
java -jar /ohta/apps/GenomeAnalysisTK.jar  -T GenotypeGVCFs -R .fasta --variant .g.vcf --variant .g.vcf -o .vcf -allSites 


#set filter for SNPs
java -jar /ohta/apps/GenomeAnalysisTK.jar -T VariantFiltration -R /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa --filterExpression 'QUAL<50.0' --filterName 'lowQUAL'  --variant ${ind}.uni.vcf -o ${ind}.filt.vcf
#only retain the SNPs you are interested in
java -jar /ohta/apps/GenomeAnalysisTK.jar -T SelectVariants -R /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -V ${ind}.filt.vcf -o ${ind}.clean.vcf -env -ef  --selectTypeToExclude INDEL --restrictAllelesTo BIALLELIC
#--maxFilteredGenotypes 0 --maxNOCALLnumber 0

