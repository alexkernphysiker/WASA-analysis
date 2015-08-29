if [ `qstat|wc -l` -gt 50 ]
then
	echo "Your task queue contains too much tasks"
	exit 1
fi
for X in `seq 45873 1 46884`
  do
	echo "Analysis for run ${X} ($1 reaction)"
	if [ -f ${RUNS_DATA}/run_${X} ]
	then
		if [ -f $PWD/../Preselection/Data$1_run_${X}.root ]
		then
			echo "Has been analysed."
		else
			if [ `qstat|wc -l` -lt 50 ]	
			then
				if [ `qstat|grep ${X}_$1|wc -l` -gt 0 ]
				then
					echo "Is being analysed."
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
					echo "HAS BEEN STARTED."
					sleep 2
				fi 
			else
				echo "Has not been started (too many jobs are running)"
			fi
		fi
	else
		echo "Has not been started (no such run)"
	fi
done

