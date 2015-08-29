if [ `qstat|grep $1|wc -l` -gt 49 ]
then
	echo "Too many jobs are running"
	exit 1
fi
for X in `seq 45873 1 46884`
  do
	if [ -f ${RUNS_DATA}/run_${X} ]
	then
		if [ -f $PWD/../Preselection/Data$1_run_${X}.root ]
		then
		else
			if [ `qstat|grep $1|wc -l` -lt 50 ]	
			then
				if [ `qstat|grep ${X}_$1|wc -l` -gt 0 ]
				then
				else
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
					echo "Analysis for run ${X} ($1 reaction) has been started"
					sleep 2
				fi 
			else
			fi
		fi
	else
	fi
done

