locList=$1
outDir=$2
outgroupName=$3


awk 'split($1,a,"."){print a[1]}' ${locList} >${outgroupName}/loc.list
mkdir ${outDir}/${outgroupName}
while read loc
do
cat ${outgroupName}/${loc}.fasta ${outgroupName}/${loc}.fasta >>${outDir}/${outgroupName}/${loc}.fasta
done < ${outgroupName}/loc.list
