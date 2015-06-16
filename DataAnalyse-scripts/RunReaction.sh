for X in ${RUNS_SEQ}
  do
     scriptname="run_${X}.sh"
     rm -f ${scriptname}
     echo "#!/bin/bash" >> ${scriptname}
     echo "#PBS -N DATA_${X}_$1" >> ${scriptname}
     echo "#PBS -l walltime=48:00:00" >> ${scriptname}

     echo >> ${scriptname}
     echo "cd $PWD/../DataAnalyse" >> ${scriptname}
     echo "./main Data$1 -fin cluster:${RUNS_DATA}/run_${X} -r ${X} -n Data$1_run_${X} -lf run_${X}.log -abort " >> ${scriptname}
     echo >> ${scriptname}
     echo "rm -f ${scriptname}" >> ${scriptname}
     chmod u+x ${scriptname}
     qsub ${scriptname}
     echo "${scriptname} generated and executed"
     sleep 2 
done
scriptname="run_mc.sh"
rm -f ${scriptname}
echo "#!/bin/bash" >> ${scriptname}
echo "#PBS -N MC__$1" >> ${scriptname}
echo "#PBS -l walltime=48:00:00" >> ${scriptname}
echo >> ${scriptname}
echo "cd $PWD/../DataAnalyse" >> ${scriptname}
echo "./main MC$1 -mode mc -fin file:$WMC_DATA/$1.wmc.data -n MC$1 -abort" >> ${scriptname}
echo >> ${scriptname}
echo "rm -f ${scriptname}" >> ${scriptname}
chmod u+x ${scriptname}
echo "${scriptname} generated"
scriptname2="run_mc_prelude.sh"
rm -f ${scriptname2}
echo "#!/bin/bash" >> ${scriptname2}
echo "#PBS -N MC_$1" >> ${scriptname2}
echo "#PBS -l walltime=48:00:00" >> ${scriptname2}
echo >> ${scriptname2} 
echo "if [-f $WMC_DATA/$1.wmc.data];" >> ${scriptname2}
echo "then" >> ${scriptname2}
echo "./$scriptname" >> ${scriptname2}
echo "else" >> ${scriptname2}
echo "cd $PWD/WMC" >> ${scriptname2}
echo "./wmc.sh $PLUTO_OUTPUT/$1.root $WMC_DATA/$1.wmc.data" >> ${scriptname2}
echo "qsub $scriptname" >> ${scriptname2}
echo "fi" >> ${scriptname2}
echo >> ${scriptname2}
echo "rm -f ${scriptname2}" >> ${scriptname2}
chmod u+x ${scriptname2}
qsub ${scriptname2}
echo "${scriptname2} generated and executed"
