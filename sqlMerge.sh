
$1 = inds
$2 = base
$3 = out

awk 'split($1,a,"."){print a[1] "" a[2]}' ${inds} >rna.sql.list
while read ind
do
echo "create table ${ind} (Locus VARCHAR(150), mapped INT);" | mysql -u felix.beaudry -D ${out}
#is file name now intosql.py?
python /ohta/felix.beaudry/scripts/reads2poly/intosql.py -i ${ind}.reads.txt -c ${ind} | mysql -u felix.beaudry -D ${out}
done < rna.sql.list

sqlMerge.py -i rna.sql.list -l ${base} -t ${out} >sqlMerging.sh
bash -i sqlMerging.sh

echo "select * from ${out};" | mysql ${out} >${out}.tsv


#sql codes
#CREATE DATABASE rhastatulus;