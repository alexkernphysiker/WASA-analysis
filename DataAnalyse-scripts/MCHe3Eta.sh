#!/bin/bash
#PBS -N analysis_task
#PBS -l walltime=12:00:00
cd ../DataAnalyse
./main MCHe3Eta Raw -mode mc -fin file:$WMC_DATA/etap.ems -n MCHe3Eta -abort 
cd ../DataAnalyse-scripts
