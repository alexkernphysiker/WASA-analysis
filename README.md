WASA Analysis
==============


Folders in this repository
==========================


General - files that are included from different projects here.

RunPluto - here is a cmake project for the application that runs pluto for simulating reactions analysed here.

DataAnalyse - make project of application that runs data preselection and Monte Carlo data analysis for reactions analysed here.

DataAnalyse-scripts - bash scripts for running preselection.

PostAnalysis - cmake project for applications that run analysis of data after preselection.



Needed environment variables
============================

ROOTSYS - path where ROOT is installed

PLUTOSYS - path where pluto is installed

WASA_ROOT - path where WMC is installed

ROOTSORTERSYS - path where RootSorter is installed

RUNS_DATA - path where data from experiment are located

WMC_DATA - path where data from WMC are located. These data are analysed further by preselection algorithm.

RUNS_SEQ - sequence for run numbers analysed by preselection algorithm. Example: `seq 46219 1 46240`

