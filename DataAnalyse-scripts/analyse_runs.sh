indir="/data0/AllRuns"
for X in ${RUNS_SEQ}
  do
     scriptname="run_${X}.sh"
     rm -f ${scriptname}
     echo "#!/bin/bash" >> ${scriptname}
     echo "#PBS -N ANALYSE_DATA_${X}" >> ${scriptname}
     echo "#PBS -l walltime=48:00:00" >> ${scriptname}

     echo >> ${scriptname}
     echo "cd $PWD/../DataAnalyse" >> ${scriptname}
     echo " ./main $1 -fin cluster:${indir}/run_${X} -r ${X} -n $1_run_${X} -lf run_${X}.log -abort " >> ${scriptname}
     echo >> ${scriptname}
     echo "rm -f ${scriptname}" >> ${scriptname}
     chmod u+x ${scriptname}
     qsub ${scriptname}
     echo "${scriptname} generated and executed"
     sleep 2 
done
