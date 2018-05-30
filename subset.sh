for set in rna rna_felix rna_josh rna_josh_males rna_josh_fem
do
mkdir subsets/$set
for chrom in auto hemi XY
do
bash /ohta/felix.beaudry/scripts/reads2poly/fasta_maker.sh ${set}.inds ${chrom}.list subsets/$set/$chrom 0.6
done
done