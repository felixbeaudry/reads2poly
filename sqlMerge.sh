
#usage bash -i sqlMerge.sh rna.inds RhastatulusBase RhastatulusExp

$1 = inds
$2 = base
$3 = out


#cat A.list H.list N.list X.list >allLoci.list
echo "create table ${base} (Loci VARCHAR(150));" | mysql -u felix.beaudry -D felix_rhastatulus
python /ohta/felix.beaudry/scripts/reads2poly/intosql_base.py -i allLoci.list -c RhastatulusBase | mysql -u felix.beaudry -D felix_rhastatulus

awk 'split($1,a,"."){print a[1] "" a[2]}' ${inds} >rna.sql.list
#awk 'split($1,a,"."){print a[1] "" a[2]}' rna.inds >rna.sql.list 

while read old
do
new=$(echo "${old}" | awk 'split($1,a,"."){print a[1] "" a[2]}')
echo "Counting reads for ${new}"
samtools index ${old}_ddpl.bam 
samtools idxstats ${old}_ddpl.bam >${new}.reads.txt
done < rna.inds

while read ind
do
echo "DROP TABLE ${ind};" | mysql -u felix.beaudry -D felix_rha
echo "CREATE TABLE ${ind} (Locus VARCHAR(150), mapped INT);" | mysql -u felix.beaudry -D felix_rha
python /ohta/felix.beaudry/scripts/reads2poly/intosql.py -i ${ind}.reads.txt -c ${ind} | mysql -u felix.beaudry -D felix_rha
done < rna.sql.list

sqlMerge.py -i rna.sql.list -l ${base} -t ${out} >sqlMerging.sh
bash -i sqlMerging.sh

echo "select * from ${out};" | mysql ${out} >${out}.tsv


#sql codes
#mysql
#CREATE DATABASE felix_rhastatulus;
#USE felix_rhastatulus;

#cat A.list H.list N.list X.list >allLoci.list