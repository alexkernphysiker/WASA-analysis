WASA Analysis (searching eta-mesic He3)
=======================================
Sources of my software used for analysis of data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.

This repository contains ubmodules with libraries from othe my repositories. You can read also their documentation.



Required software
=================

ROOT - framework for calulations

Pluto - library for Monte Carlo simulation of reaction products. Uses ROOT.

WMC - software for Monte Carlo simulation of registration of previously simulated particles by WASA detector. Uses ROOT and Pluto.

RootSorter - framework for data preselection (both WMC-generated and obtained from the measurements). Uses ROOT.

gnuplot - software for plotting. Is used by some applications performing final analysis.



Needed environment variables
============================

Needed software:

ROOTSYS - path where ROOT is installed

PLUTOSYS - path where pluto is installed

WASA_ROOT - path where WMC is installed

ROOTSORTERSYS - path where RootSorter is installed

This project:

PLUTO_OUTPUT - path where pluto files are stored

RUNS_DATA - path where data from experiment are located

WMC_DATA - path where data from WMC are located. These data are analysed further by preselection algorithm.



Directories
===========

Preselection - Makefile project of software that provides data preselection. 
Is compiled and runs on wasa00 server. 
ATTENTION! Includes config.h from upper directory and contains relative symlinks to other used files from this repo.

DataFiles - additional information required for data preselection.

wasa_scripts - directory with scripts and additional data for configuring used WASA software and running WASA Monte Carlo and data preselection.

lib - contains one used submodule (FitGen) and sources of some algoritms implemented as libraries

test - will contain unit tests for parts of this project :)

./ - contains main cmake project with applications for: running pluto simulations of needed reactions, providing reconstruction fits and analysis of data obtained after preselection



Compiling and running preselection (on wasa00)
============================================

git clone https://github.com/alexkernphysiker/WASA-analysis.git

cd WASA-analysis

git submodule update --init --recursive

cd wasa_scripts/config

./config.sh #configure WMC and rootsorter for current experimental setup

cd ../../Preselection

make # compile preselection software

cd ../wasa_scripts

Here you can run scripts for starting either WMC or data preselection



Compiling other programs (locally)
========================

git clone https://github.com/alexkernphysiker/WASA-analysis.git

cd WASA-analysis

git submodule update --init --recursive

cmake .

make

Here you can run application for analysis

