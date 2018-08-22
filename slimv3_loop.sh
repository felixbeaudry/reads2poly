for typ in n p
do
        for sel in {0..300..5}
        do
                for pro in {0..10}
                do
                        python /ohta/felix.beaudry/scripts/reads2poly/slim_param_v3.py -s ${sel} -p ${pro} -t ${typ} >slims/${typ}_${sel}_${pro}.slim
                        echo "mean se slope" >${typ}_${sel}_${pro}.txt
                        for sim in {0..100}
                        do
                        echo "${typ} ${sel} ${pro} ${sim}"
                        /ohta/apps/SLiM/bin/slim slims/${typ}_${sel}_${pro}.slim | awk '$1 ~ "ratio"{print $2}' >sims/${typ}_${sel}_${pro}_${sim}.txt
                        python /ohta/felix.beaudry/scripts/reads2poly/slim_stat.py -i sims/${typ}_${sel}_${pro}_${sim}.txt >>${typ}_${sel}_${pro}.txt
                        done
                done
        done
done