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
    const auto runs=PresentRuns("L");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"ppn");
    Plot<double>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_0"),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_0"),"ppn_{sp}")<<"set key on"
    <<"set title 'Planarity. MC'"<<"set yrange [0:]";
    Plot<double>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_1"),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_1"),"ppn_{sp}")<<"set key on"
    <<"set title 'Planarity. MC. Cut'"<<"set yrange [0:]";
    Plot<double>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"pair_phi_diff_0"))
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"pair_phi_diff_1"))
    <<"set title 'Planarity. Data "+runmsg+"'"<<"set yrange [0:]";

    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. MC pd'";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. MC ppn_{sp}'";
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. Data "+runmsg+"'";


    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_21"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. MC pd'";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_21"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. MC ppn_{sp}'";
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_21"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. Data "+runmsg+"'";
    Plot<double>()
    .Line(Hist(MC,"pd",{"Histograms","elastic"},"theta_sum_21-AllBins").toLine(),"pd")
    .Line(Hist(MC,"ppn_qf",{"Histograms","elastic"},"theta_sum_21-AllBins").toLine(),"ppn_{sp}")
    <<"set title 'MC'"<<"set key on"<<"set yrange [0:]"<<"set xrange [0:210]";
    Plot<double>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"theta_sum_21-AllBins"))
    <<"set title 'Data "+runmsg+"'"<<"set key on"<<"set yrange [0:]"<<"set xrange [0:210]";


    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. MC pd'";
    PlotHist2d<double>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. MC ppn_{sp}'";
    PlotHist2d<double>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Theta vs Theta. Data "+runmsg+"'";
    Plot<double>()
    .Line(Hist(MC,"pd",{"Histograms","elastic"},"theta_sum_22-AllBins").toLine(),"pd")
    .Line(Hist(MC,"ppn_qf",{"Histograms","elastic"},"theta_sum_22-AllBins").toLine(),"ppn_{sp}")
    <<"set title 'MC'"<<"set key on"<<"set yrange [0:]"<<"set xrange [0:210]";
    Plot<double>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"theta_sum_22-AllBins"))
    <<"set title 'Data "+runmsg+"'"<<"set key on"<<"set yrange [0:]"<<"set xrange [0:210]";
}
