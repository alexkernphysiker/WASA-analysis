// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
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
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "he3-forward-kinematics");
    const string He3eta_msg = "He3eta-gg";
    const string bound_msg = "bound1-2g";
    auto QBins = Hist(MC, He3eta_msg, {"Histograms", "He3Forward_Reconstruction"}, "0-Reference");
    for (size_t bin_num = 0, bin_count = QBins.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = QBins[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream() << "Q in [" << setprecision(3)
                            << Q.min() << "; " << Q.max() << "] MeV").str();
        auto kin_v = Hist2d(MC, He3eta_msg,
        {"Histograms", "He3Forward_Vertices"},
        string("Kinematic-vertex-Bin-") + to_string(bin_num)).Scale(5, 5);
        PlotHist2d(sp2,"",5).Distr(kin_v, Qmsg) << "set key on" << "set xlabel 'E_k, GeV'"
                                             << "set ylabel 'theta, deg'" << "set xrange [0.2:0.4]";
        auto kin_bound = Hist2d(MC, bound_msg,
        {"Histograms", "He3Forward_Reconstruction"},
        string("Kinematic-reconstructed-Bin-") + to_string(bin_num)
                            ).Scale(5, 5);
        PlotHist2d(sp2, Q.Contains(21) ? "He3forward-kin-bound-above" :(Q.Contains(-10) ? "He3forward-kin-bound-below":""),5).Distr(kin_bound)
                << "set xlabel 'E_k, GeV'"
                << "set ylabel 'theta, deg'"
                << "set xrange [0.2:0.4]"
                << "set title 'Monte Caro, " + Qmsg + "'";
        auto kin_mc = Hist2d(MC, He3eta_msg,
        {"Histograms", "He3Forward_Reconstruction"},
        string("Kinematic-reconstructed-Bin-") + to_string(bin_num)
                            ).Scale(5, 5);
        PlotHist2d(sp2, Q.Contains(21) ? "He3forward-kin-mc" : "",5).Distr(kin_mc)
                << "set xlabel 'E_k, GeV'"
                << "set ylabel 'theta, deg'"
                << "set xrange [0.2:0.4]"
                << "set title 'Monte Caro, " + Qmsg + "'";
        auto data_hist = Hist2d(DATA, "All",
        {"Histograms", "He3Forward_Reconstruction"},
        string("Kinematic-reconstructed-Bin-") + to_string(bin_num)
                               ).Scale(5, 5);
        PlotHist2d(sp2, Q.Contains(21) ? "He3forward-kin-data" : "",5).Distr(data_hist)
                << "set xlabel 'E_k, GeV'"
                << "set ylabel 'theta, deg'"
                << "set xrange [0.2:0.4]"
                << "set title 'Data (with correction), " + Qmsg + "'";
        const auto xcut = (kin_mc.X().size() * 2) / 5;
        const auto &xC = kin_mc.X()[xcut].min();
        const hist<> ymc = kin_mc.CutY(xcut), ydata = data_hist.CutY(xcut);
        Plot("",7).Hist(ymc / ymc.TotalSum().val(), "MC").Hist(ydata / ydata.TotalSum().val(), "DATA")
                << "set key on" << "set yrange[0:0.1]" << "set ylabel ''"<<"set title'"+Qmsg+"'"
                << "set xlabel 'theta,deg (E=" + to_string(xC) + " GeV)'"<<"set xrange[4:9]";
        const auto ycut = kin_mc.Y().size() / 3;
        const auto &yC = kin_mc.Y()[ycut].min();
        const hist<> xmc = kin_mc.CutX(ycut), xdata = data_hist.CutX(ycut);
        Plot("",7).Hist(xmc / xmc.TotalSum().val(), "MC").Hist(xdata / xdata.TotalSum().val(), "DATA")
                << "set key on" << "set yrange[0:]" << "set ylabel ''"<<"set title'"+Qmsg+"'"
                << "set xlabel 'E_k, GeV (theta=" + to_string(yC) + " deg)'"
                << "set xrange [0.2:0.4]";
        const auto ycut2 = kin_mc.Y().size() * 2 / 3;
        const hist<> xmc2 = kin_mc.CutX(ycut2), xdata2 = data_hist.CutX(ycut2);
        Plot("",7).Hist(xmc2 / xmc2.TotalSum().val(), "MC").Hist(xdata2 / xdata2.TotalSum().val(), "DATA")
                << "set key on" << "set yrange[0:]" << "set ylabel ''"<<"set title'"+Qmsg+"'"
                << "set xlabel 'E_k, GeV (theta=" + to_string(yC) + " deg)'"
                << "set xrange [0.2:0.4]";
    }
    return 0;
}
