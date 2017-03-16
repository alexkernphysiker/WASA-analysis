// this file is distributed under 
// GPL license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/functions.h>
#include <math_h/sigma.h>
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
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"beam_momenta");
	string He3eta_msg="He3eta";
	auto QBins=Hist(MC,He3eta_msg,{"Histograms","He3Forward_Reconstruction"},"0-Reference");
	for(size_t bin_num=0,bin_count=QBins.size();bin_num<bin_count;bin_num++){
		string Qmsg="Q in ["+to_string(QBins[bin_num].X().min())+":"+to_string(QBins[bin_num].X().max())+"] MeV";
		auto kin_v=Hist2d(MC,He3eta_msg,{"Histograms","He3Forward_Vertices"},string("Kinematic-vertex-Bin-")+to_string(bin_num)).Scale(5,5);
		PlotHist2d<double>(sp2).Distr(kin_v,Qmsg)<<"set key on" << "set xlabel 'E_k, GeV'" 
		<< "set ylabel 'theta, deg'"<< "set xrange [0.2:0.4]";
		auto kin_mc=Hist2d(MC,He3eta_msg,{"Histograms","He3Forward_Reconstruction"},string("Kinematic-reconstructed-Bin-")+to_string(bin_num)).Scale(5,5);
		PlotHist2d<double>(sp2).Distr(kin_mc,Qmsg)<<"set key on" <<"set xlabel 'E_k, GeV'"
		<<"set ylabel 'theta, deg'"<< "set xrange [0.2:0.4]";
		auto data_hist=Hist2d(DATA,"",{"Histograms","He3Forward_Reconstruction"},string("Kinematic-reconstructed-Bin-")+to_string(bin_num)).Scale(5,5);
		PlotHist2d<double>(sp2).Distr(data_hist,Qmsg)<<"set key on" <<"set xlabel 'E_k, GeV'"
		<<"set ylabel 'theta, deg'"<< "set xrange [0.2:0.4]";
		const auto xcut=(kin_mc.X().size()*2)/5;
		const auto&xC=kin_mc.X()[xcut].val();
		const hist<double> ymc=kin_mc.CutY(xcut),ydata=data_hist.CutY(xcut);
		Plot<double>().Hist(ymc/ymc.TotalSum().val(),"MC").Hist(ydata/ydata.TotalSum().val(),"DATA")
		<<"set key on"<<"set yrange[0:]"<<"set ylabel ''"
		<<"set xlabel 'theta,deg (E="+to_string(xC)+" GeV)'";
		const auto ycut=kin_mc.Y().size()/3;
		const auto&yC=kin_mc.Y()[ycut].val();
		const hist<double> xmc=kin_mc.CutX(ycut),xdata=data_hist.CutX(ycut);
		Plot<double>().Hist(xmc/xmc.TotalSum().val(),"MC").Hist(xdata/xdata.TotalSum().val(),"DATA")
		<<"set key on"<<"set yrange[0:]"<<"set ylabel ''"
		<<"set xlabel 'E_k, GeV (theta="+to_string(yC)+" deg)'"
		<< "set xrange [0.2:0.4]";
	}
	return 0;
}
