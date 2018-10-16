##sqlMerge.sh Felix Beaudry Oct 10th 2018
##usage: bash -i sqlMerge.sh rna.inds RhastatulusBase RhastatulusExp

$1 = inds
$2 = base
$3 = out

##Make a table with all the loci onto which to left join
#cat A.list H.list N.list X.list >allLoci.list
echo "create table ${base} (Loci VARCHAR(150));" | mysql -u felix.beaudry -D felix_rhastatulus
python /ohta/felix.beaudry/scripts/reads2poly/intosql_base.py -i allLoci.list -c RhastatulusBase | mysql -u felix.beaudry -D felix_rhastatulus

##remove periods from individual names with periods, otherwise SQL will not accept
awk 'split($1,a,"."){print a[1] "" a[2]}' ${inds} >rna.sql.list

##Count reads and remove periods from individual's names
while read old
do
new=$(echo "${old}" | awk 'split($1,a,"."){print a[1] "" a[2]}')
echo "Counting reads for ${new}"
samtools index ${old}_ddpl.bam 
samtools idxstats ${old}_ddpl.bam >${new}.reads.txt
done < rna.inds

##Insert read counts into SQL
while read ind
do
##echo "DROP TABLE ${ind};" | mysql -u felix.beaudry -D felix_rha
echo "CREATE TABLE ${ind} (Locus VARCHAR(150), mapped INT);" | mysql -u felix.beaudry -D felix_rha
python /ohta/felix.beaudry/scripts/reads2poly/intosql.py -i ${ind}.reads.txt -c ${ind} | mysql -u felix.beaudry -D felix_rha
done < rna.sql.list

##use list of names to make a SQL merge command
python /ohta/felix.beaudry/scripts/reads2poly/sqlMerge.py -i rna.sql.list -l ${base} -t ${out}
##echo "DROP TABLE ${ind};" | mysql -u felix.beaudry -D felix_rha
bash -i sqlMerging.sh

##export merged database
echo "select * from ${out};" | mysql -u felix.beaudry -D felix_rha >${out}.tsv


####SQL cheatsheet####
#CREATE DATABASE felix_rhastatulus;
#USE felix_rhastatulus;
