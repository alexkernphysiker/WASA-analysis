#/bin/bash
cp FDEdep2Ekin_3He $ROOTSORTERSYS/wasa/Database
cp Calibration/* $ROOTSORTERSYS/wasa/Database/PreviousRuns/
mkdir $WASA_ROOT/alig/m4/B009
cp WMC_Geometry/* $WASA_ROOT/alig/m4/B009
cp al4cosy0_mc_009.dat042.m4 $WASA_ROOT/alig/m4
cd $WASA_ROOT/alig/m4 
mv B009/mkalig.sh .
./mkalig.sh

