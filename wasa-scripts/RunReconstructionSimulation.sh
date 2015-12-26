scriptname="run_mc.sh"
rm -f ${scriptname}
echo "#!/bin/bash" >> ${scriptname}
echo "#PBS -N RE_$1" >> ${scriptname}
echo "#PBS -l walltime=48:00:00" >> ${scriptname}
echo >> ${scriptname}
echo "cd $PWD/../Preselection" >> ${scriptname}
echo "./main RE_$1 -mode mc -fin file:$WMC_DATA/$1.wmc.data -n RE$1 -abort" >> ${scriptname}
echo >> ${scriptname}
echo "rm -f $PWD/${scriptname}" >> ${scriptname}
chmod u+x ${scriptname}
qsub ${scriptname}
echo "${scriptname} generated and executed"
qstat

