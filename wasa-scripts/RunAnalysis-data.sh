if [ `qstat|grep $1|wc -l` -gt 0 ]
then
	echo "The task queue allready contains tasks for this reaction"
	exit 1
fi
for X in `seq 45873 1 46884`
  do
	if [ -f ${RUNS_DATA}/run_${X} ]
	then
		if [ -f $PWD/../Preselection/Data$1_run_${X}.root ]
		then
	echo "Data$1_run_${X}.root already exists"
		else
			if [ `qstat|wc -l` -lt 50 ]	
			then
     scriptname="run_${X}.sh"
     rm -f ${scriptname}
     echo "#!/bin/bash" >> ${scriptname}
     echo "#PBS -N DATA_${X}_$1" >> ${scriptname}
     echo "#PBS -l walltime=48:00:00" >> ${scriptname}

     echo >> ${scriptname}
     echo "cd $PWD/../Preselection" >> ${scriptname}
     echo "./main Data_$1 -fin cluster:${RUNS_DATA}/run_${X} -r ${X} -n Data$1_run_${X} -lf run_${X}.log -abort " >> ${scriptname}
     echo >> ${scriptname}
     echo "rm -f $PWD/${scriptname}" >> ${scriptname}
     chmod u+x ${scriptname}
     qsub ${scriptname}
     echo "${scriptname} generated and executed"
     sleep 2 
			else
			echo "Too many jobs"
			fi
		fi
	else
	echo "no file ${RUNS_DATA}/run_${X}"
	fi
done

