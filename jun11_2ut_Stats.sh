for pop in FL #NCsub #FLNC
do
for chrom in Y #X 
do
cp ${pop}male${chrom}.pop subset7/loose/NXY/rothschildianus/${pop}/${pop}male${chrom}.pop
perl /ohta/felix.beaudry/scripts/reads2poly/Polymorphurama_v3.1.pl -i subset7/loose/NXY/rothschildianus/${pop}/ -p ${pop}male${chrom} -S loose -o rothschildianus -C N${chrom}
done
done



cp TXFL.pop subsets2/XXY/TXFL.pop
perl /ohta/felix.beaudry/scripts/reads2poly/Polymorphurama_v3.1.pl -i subsets2/XXY/ -p TXFL -S loose -o NA -C Y -c Ychrom






cp X2Y.pop subset7/loose/NXY/both/X2Y.pop

perl /ohta/felix.beaudry/scripts/reads2poly/Polymorphurama_v3.1.pl -i subset7/loose/NXY/both/ -p X2Y -S loose -o bucephalophorus -C NX2Y

cp X2Y.pop subsets2/XXY/X2Y.pop
perl /ohta/felix.beaudry/scripts/reads2poly/Polymorphurama_v3.1.pl -i subsets2/XXY/ -p X2Y -S loose -o bucephalophorus -C X2Y