// this file is distributed under 
// GPL license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include "he3eta.h"
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward_cuts");
	auto MakePlots=[](const histsource HS,const string& R, const uint max_z){
		string title="Data";
		if(HS==MC){
			title="Simulation "+R;
			Plot<double>()
			.Line(Hist(HS,R,{"Histograms","He3Forward_Reconstruction"},"0-Reference").toLine(),"Simulated events")
			.Line(Hist(HS,R,{"Histograms","He3Forward_Reconstruction"},"2-ThetaIsAccepted").toLine(),"Reconstructable forward tracks")
			.Line(Hist(HS,R,{"Histograms","He3Forward_Reconstruction"},"4-GeomCut").toLine(),"identified as ^3He")
			<< "set key on"<< "set title '"+title+"'"
			<< "set yrange [0:4000000]"
			<< "set xlabel 'Q, MeV'"
			<< "set ylabel 'Events count'";
		}
		auto SP2=[max_z](){
			return PlotHist2d<double>(sp2)
			<<"set key on"
			<<"set xrange [0:0.3]"<<"set xlabel 'E_{FRH1}, GeV'"
			<<"set yrange [0:0.03]"<<"set ylabel 'E_{FTH1}, GeV'"
			<<"set zrange [0:"+to_string(max_z)+"]";
		};
		auto THD=[](){
			return Plot<double>()
			<<"set key on" 
			<< "set xlabel 'Phi, deg'"<<"set xrange [-180:180]"
			<< "set ylabel 'events, count'"<< "set yrange [0:]";
		};
		string fpc="set title 'Reconstructable tracks. "+title+"'";
		SP2().Distr(Hist2d(HS,R,{"Histograms","He3Forward_Reconstruction"},string("2-ThetaIsAccepted-FTH1-vs-FRH1")).Scale(3,3))<<fpc;
		THD().Hist(Hist(HS,R,{"Histograms","He3Forward_Debug"},"2-PhiDistribution-AllBins").Scale(30))<<fpc;
		string cut="set title 'Identified as 3He. "+title+"'";
		SP2().Distr(Hist2d(HS,R,{"Histograms","He3Forward_Reconstruction"},string("4-GeomCut-FTH1-vs-FRH1")).Scale(3,3))<<cut;
		THD().Hist(Hist(HS,R,{"Histograms","He3Forward_Debug"},"4-PhiDistribution-AllBins").Scale(30))<<cut;
	};
	MakePlots(MC,"He3eta",500000);
	MakePlots(MC,"He3pi0pi0",500000);
	MakePlots(MC,"He3pi0pi0pi0",500000);
	MakePlots(DATA,"L",50000);
	MakePlots(MC,"He3pi0",500000);
}
