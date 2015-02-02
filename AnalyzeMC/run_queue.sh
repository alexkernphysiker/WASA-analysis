#!/bin/bash
#PBS -N mc_dd_3Hen
#PBS -l walltime=12:00:00
./main -mode mc -fin $1 -n analysisjob -abort 

