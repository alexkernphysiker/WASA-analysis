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
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"Compl_mag_0"),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"Compl_mag_0"),"ppn_qf")<<"set key on";
    Plot<double>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"Compl_mag_2"),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"Compl_mag_2"),"ppn_qf")<<"set key on";
    Plot<double>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"Compl_mag_0"))
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"Compl_mag_1"))
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"Compl_mag_2"));
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:800]";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_2"))
    <<"set zrange [0:]";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_2"))
    <<"set zrange [0:800]";
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:800]";
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_2"))
    <<"set zrange [0:800]";
    Plot<double>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"theta_sum_1"))
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"theta_sum_2"));
}
