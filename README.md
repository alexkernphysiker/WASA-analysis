WASA Analysis
=============
Sources of my software used for analysis of data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.
This repository contains code that does not require RootSorter.

All files are distributed under GPL license


Required software
=================

Framework for calulations that is required to read the input data

	ROOT

Software for plotting. Is used by some applications performing final analysis.

	gnuplot

I compile this code with the version of gcc new enough to support C++17 standard.

Needed environment variables
============================

Path containing raw-data analysis results

	WASA_OUTPUT_DATA

Path to store final results of analysis

	POST_ANALYSIS_DATA


Compiling
=========

	git clone https://github.com/alexkernphysiker/WASA-analysis.git
	cd WASA-analysis
	git submodule update --init --recursive
	cmake .
	make
