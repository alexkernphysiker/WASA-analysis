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
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_0"),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_0"),"ppn_qf")<<"set key on";
    Plot<double>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_1"),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_1"),"ppn_qf")<<"set key on";
    Plot<double>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"pair_phi_diff_0"))
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"pair_phi_diff_1"));
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_21"))
    <<"set zrange [0:]";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_21"))
    <<"set zrange [0:]";
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]";
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_21"))
    <<"set zrange [0:]";
    Plot<double>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"theta_sum_21"))
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"theta_sum_21"));
    Plot<double>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"theta_sum_21"));
}
