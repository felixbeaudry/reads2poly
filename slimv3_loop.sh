for typ in n p
do
	for sel in {0..500}
	do
		for pro in {0..10}
		do
			python /ohta/felix.beaudry/scripts/reads2poly/slim_param_v3.py -s 0.${sel} -p 0.${pro} -t ${typ} >${typ}_${sel}_${pro}.slim 
			echo "mean se slope" >${typ}_${sel}_${pro}.txt
			for sim in {0..100}
			do
			echo "${typ} ${sel} ${pro} ${sim}"
			/ohta/apps/SLiM/bin/slim ${typ}_${sel}_${pro}.slim | awk '$1 ~ "ratio"{print $2}' >${typ}_${sel}_${pro}_${sim}.txt 
			#perl /ohta/felix.beaudry/scripts/reads2poly/slim_stat.pl ${typ}_${sel}_${pro}_${sim}.txt >>${typ}_${sel}_${pro}.txt
			#rm ${typ}_${sel}_${pro}_${sim}.txt
			done
		done
	done
done

