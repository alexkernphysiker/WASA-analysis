scriptname="run_wmc.sh"
rm -f ${scriptname}
echo "#!/bin/bash" >> ${scriptname}
echo "#PBS -N WMC_$1" >> ${scriptname}
echo "#PBS -l walltime=48:00:00" >> ${scriptname}
echo >> ${scriptname}
echo "cd $PWD/WMC" >> ${scriptname}
echo "./wmc.sh $PLUTO_OUTPUT/$1.root $WMC_DATA/$1.wmc.data" >> ${scriptname}
echo >> ${scriptname}
echo "rm -f $PWD/${scriptname}" >> ${scriptname}
chmod u+x ${scriptname}
qsub ${scriptname}
echo "${scriptname} generated and executed"
