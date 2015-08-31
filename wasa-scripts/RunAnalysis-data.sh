for X in `seq 45873 1 46884`; do
	if [ -f ${RUNS_DATA}/run_${X} ]; then
		if [ ! -f $PWD/../Preselection/Data$1_run_${X}.root ]; then
			if [ `qstat|grep $1|grep "R \|Q "|wc -l` -lt 50 ]; then
				if [ `qstat|grep ${X}_$1|wc -l` -lt 1 ]; then
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
				echo "We have already enough jobs running for this reaction"
				exit 0
			fi
		fi
	fi
done
qstat

