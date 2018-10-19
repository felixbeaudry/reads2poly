/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/blastn -task megablast -query /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -db /ohta/felix.beaudry/assemblies/bucephalophorusTranscriptome/RB1_transcriptome_ref.fa -num_threads 10 -max_target_seqs 1 -outfmt 6 -dust no -gapopen 0 -gapextend 0 >bucephalophorus.blast 2>bucephalophorus.blast.err

/ohta/aplatts/data/apps/ncbi-blast-2.6.0+/bin/blastn -task megablast -query /ohta/felix.beaudry/assemblies/hastTranscriptome/NCF1_combined_ref_CDS_noStop.fa -db /ohta/felix.beaudry/assemblies/rothTranscriptome/roth_female_ref_genewise.fa  -num_threads 10 -max_target_seqs 1 -outfmt 6 -dust no -gapopen 0 -gapextend 0 >rothschildianus.blast 2>rothschildianus.blast.err

awk '{print $1"\t"$2}' bucephalophorus.blast | uniq >bucephalophorus.orthologs
awk '{print $1"\t"$2}' rothschildianus.blast | uniq >rothschildianus.orthologs