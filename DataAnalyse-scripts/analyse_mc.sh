#!/bin/bash
#PBS -N ANALYSE_WMC
#PBS -l walltime=12:00:00
cd ~/WASA-analysis/DataAnalyse
./main $1 -mode mc -fin file:$WMC_DATA/etap.ems -n $1 -abort 
cd ../DataAnalyse-scripts
