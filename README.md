WASA Analysis
=============
Sources of my software used for analysis of data obtained from the experiment WASA-at-COSY on searching eta-mesic 3He in May 2014.


Required software
=================
	ROOT 
framework for calulations

	gnuplot
software for plotting. Is used by some applications performing final analysis.


Needed environment variables
============================
	ROOTSYS 
path where ROOT is installed from WASA-libs/config.h


Compiling (locally)
===================

	git clone https://github.com/alexkernphysiker/WASA-analysis.git
	cd WASA-analysis
	git submodule update --init --recursive
	cmake .
	make
