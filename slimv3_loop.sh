for typ in n p
do
	for sel in {0..500}
	do
		for pro in {0..10}
		do
			python /ohta/felix.beaudry/scripts/reads2poly/slim_param_v3.py -s 0.${sel} -p 0.${pro} -t ${typ} >${typ}_${sel}_${pro}.slim 
			for sim in {0..100}
			do
			/ohta/apps/SLiM/bin/slim ${typ}_${sel}_${pro}.slim | awk '$1 ~ "ratio"{print $2}' >${typ}_${sel}_${pro}_${sim}.txt
			# slim | awk | print mean se slope for each sim > .txt
			done
		done
	done
done

