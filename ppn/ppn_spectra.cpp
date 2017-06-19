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
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"ppn");
    Plot<double>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"Compl_mag_before"))
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"Compl_mag_after"));
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t"));
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_qf"));
}
