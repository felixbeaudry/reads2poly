



for outgroup in rothschildianus bucephalophorus
do

for chrom in X Y 
do
mkdir subset4/${outgroup}/${chrom}

while read loc
do

#rm subset4/${outgroup}/${chrom}/${loc}.fasta
samtools faidx subsets2/phase/${loc}.fasta ${outgroup} >>subset4/${outgroup}/${chrom}/${loc}.fasta

while read ind
do
samtools faidx subsets2/phase/${loc}.fasta ${ind}_${chrom}chrom >>subset4/${outgroup}/${chrom}/${loc}.fasta
done <male.inds

done <XY.list

done

done


for outgroup in rothschildianus bucephalophorus
do

for chrom in X Y 
do

cp pop.pop subset4/${outgroup}/${chrom}/pop.pop

perl /ohta/felix.beaudry/scripts/reads2poly/Polymorphurama_v3.pl -i subset4/${outgroup}/${chrom}/ -p pop -S sub4 -o ${outgroup} -C ${chrom} -c ${chrom}chrom

done

done

##

for outgroup in rothschildianus bucephalophorus
do

for chrom in A H 
do
mkdir subset4/${outgroup}/${chrom}

while read loc
do

rm subset4/${outgroup}/${chrom}/${loc}.fasta
samtools faidx subsets2/${chrom}/${loc}.fasta ${outgroup} >>subset4/${outgroup}/${chrom}/${loc}.fasta

while read ind
do
samtools faidx subsets2/${chrom}/${loc}.fasta ${ind}_1 >>subset4/${outgroup}/${chrom}/${loc}.fasta
samtools faidx subsets2/${chrom}/${loc}.fasta ${ind}_2 >>subset4/${outgroup}/${chrom}/${loc}.fasta
done <male.inds

done <${chrom}.list

done

done


for outgroup in rothschildianus bucephalophorus
do

for chrom in A H 
do

cp pop.pop subset4/${outgroup}/${chrom}/pop.pop

perl /ohta/felix.beaudry/scripts/reads2poly/Polymorphurama_v3.pl -i subset4/${outgroup}/${chrom}/ -p pop -S sub4 -o ${outgroup} -C ${chrom}

done

done