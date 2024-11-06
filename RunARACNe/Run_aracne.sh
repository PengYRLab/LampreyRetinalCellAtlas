
#! /bin/bash



# Please follow the tutorialat https://github.com/califano-lab/ARACNe-AP.


#-------------- run aracne threshold
command="java -Xmx1G -jar $aracne -e ${dir_mat}/${mat}.dat  -o ${mat}/${reg_type} --tfs ${reg} --pvalue 1E-8 --seed 1 --calculateThreshold"

echo $command | qsub -l h_rt=00:20:00 -l h_data=2G -l h_vmem=10G -N aracne_threshold_${mat}_$reg_type -j yes -o ${mat}/threshold_${mat}_${reg_type}.log -cwd



#-------------- check the arace_threshod

for i in $(seq 1 1000000)
do
   if (( $(qstat -u junqiang| grep "aracne_thr"| wc -l) > 0  )) 
   then
   echo " running the aracne threshold jobs:$(qstat -u junqiang|grep "aracne_thr"| wc -l )"
    echo "Job submission pasued"
   sleep 10s
   else
  # echo "aracne threshold is finished!"
   break
   fi
done

echo "arance threshold is finished"



#-------------run aracne reg


bootstrap_i=0


for i in $(seq 1 1000000)
do

   if (( $bootstrap_i < $bootstrap_n )) && (( $(qstat -u junqiang| grep "aracne_reg"| wc -l) > $job_max || $(qstat -u junqiang| grep "aracne_reg"| wc -l)== $job_max )) ; then
   echo " running the aracne reg jobs:$(qstat -u junqiang|grep "aracne_reg"| wc -l )"
   echo "Job submission pasued"
   sleep 10s


   elif (( $bootstrap_i < $bootstrap_n )) &&  (($(qstat -u junqiang| grep "aracne_reg"| wc -l) < $job_max )) ; then

#run aracne tfs

((bootstrap_i = bootstrap_i +1))

command="java -Xmx1G -jar $aracne -e ${dir_mat}/${mat}.dat  -o ${mat}/${reg_type} --tfs ${reg} --pvalue 1E-8 --threads 8 --seed ${bootstrap_i}"

echo $command | qsub -l h_rt=00:05:00 -l h_data=2G -l h_vmem=10G  -N aracne_reg_${mat}_${reg_type} -j yes -o ${mat}/${reg_type}_${bootstrap_i}.log -cwd

echo "submitting job $bootstrap_i"

   else
   break
   fi

done



#-------------- check the arace_reg

for i in $(seq 1 1000000)
do
   if (( $(qstat -u junqiang| grep "aracne_reg"| wc -l) > 0  )) 
   then
   echo " running the aracne threshold jobs:$(qstat -u junqiang|grep "aracne_thr"| wc -l )"
    echo "Job submission pasued"
   sleep 10s
   else
   break
   fi
done

echo "arance reg is finished"




command="java -Xmx1G -jar $aracne  -o ${mat}/${reg_type} --consolidate"


echo $command | qsub -l h_rt=00:10:00 -l h_data=2G -l h_vmem=10G  -N aracne_consolidate_${mat}_${reg_type} -j yes -o ${mat}/${reg_type}_consolidate.log -cwd


for i in $(seq 1 1000000)
do
   if (( $(qstat -u junqiang| grep "aracne_consolidate"| wc -l) > 0  )) 
   then
   echo " running the aracne threshold jobs:$(qstat -u junqiang|grep "aracne_thr"| wc -l )"
    echo "Job submission pasued"
   sleep 10s
   else
   break
   fi
done

echo "arance consolidate is finished"




