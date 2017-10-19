WASA Analysis
=============
Sources of my software used for analysis of data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.
All files are distributed under GPL license


Required software
=================
	ROOT 
framework for calulations that is required to read the input data

	gnuplot
software for plotting. Is used by some applications performing final analysis.


Needed environment variables
============================

    WASA_OUTPUT_DATA
path to store raw analysis results

    POST_ANALYSIS_DATA
path to store final results of analysis


Compiling
=========

	git clone https://github.com/alexkernphysiker/WASA-analysis.git
	cd WASA-analysis
	git submodule update --init --recursive
	cmake .
	make
