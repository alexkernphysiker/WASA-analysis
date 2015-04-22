
#PBS -q mc
#PBS -N WMC
#PBS -l walltime=2:00:00:00
cd WMC
./wmc.sh $PLUTO_OUTPUT/He3Eta_gg_.root $WMC_DATA/He3Eta.wmc.data
