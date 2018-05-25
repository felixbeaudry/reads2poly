

for outgroup in rothschildianus bucephalophorus
do
mkdir subsets/XYphased_${outgroup}
awk 'split($1,a,"."){print a[1]}' goodsex.fasta.list >subsets/XYphased_${outgroup}/goodsex.list
while read loc
do
cat XYphased/goodsex/${loc}.fasta ${outgroup}/${loc}.fasta >subsets/XYphased_${outgroup}/${loc}.fasta
done <subsets/XYphased_${outgroup}/goodsex.list
done

perl /ohta/felix.beaudry/scripts/reads2poly/polymorphurama_interpop.pl fasta XYphased_rothschildianus pop rothschildianus _Ychrom
perl /ohta/felix.beaudry/scripts/reads2poly/polymorphurama_interpop.pl fasta XYphased_rothschildianus pop rothschildianus _Xchrom



