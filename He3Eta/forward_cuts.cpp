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
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "he3eta-forward-cuts");
    auto MakePlots = [](const histsource HS, const string & R, const uint max_z) {
        string title = "";
        if (HS == MC) {
            title = "Simulation " + R;
            Plot("",5)
            .Hist(Hist(HS, R, {"Histograms", "He3Forward_Reconstruction"}, "0-Reference"), "Simulated events")
            .Hist(Hist(HS, R, {"Histograms", "He3Forward_Reconstruction"}, "2-ThetaIsAccepted"), "Reconstructable forward tracks")
            .Hist(Hist(HS, R, {"Histograms", "He3Forward_Reconstruction"}, "4-GeomCut"), "identified as ^3He")
                    << "set key on" << "set title '" + title + "'"
                    << "set xlabel 'Q, MeV'"
                    << "set ylabel 'Events count'"
                    << "set yrange [0:]"
                    << "set xrange [0:30]";
        }
        auto SP2 = [max_z](const string &name = "") {
            return PlotHist2d(sp2, name,5)
		    << "set xrange [0:0.3]"
		   << "set ytics ('0' 0,'10' 0.01,'20' 0.02,'30' 0.03) rotate by 90 right"
		    >> "set ytics"
		   << "set xtics ('0' 0,'100' 0.1,'200' 0.2,'300' 0.3)"
                    << "set xlabel 'E_{FRH1}, MeV'"
                   << "set yrange [0:0.03]" << "set ylabel 'E_{FTH1}, MeV'"
                   << "set zrange [0:" + to_string(max_z) + "]";
        };
        auto THD = [](const string &name = "") {
            return Plot(name,5)
                   << "set key on"
                   << "set xlabel 'Phi, deg'" << "set xrange [-180:180]"
                   << "set ylabel 'events, count'" << "set yrange [0:]";
        };
        const string fpc = "set title 'Reconstructable tracks " + title + "'";
        SP2((HS == DATA) ? "SP2-data-cut0" : "").Distr(Hist2d(HS, R, {"Histograms", "He3Forward_Reconstruction"}, string("2-ThetaIsAccepted-FTH1-vs-FRH1")).Scale(3, 3)) << fpc;
        THD().Hist(Hist(HS, R, {"Histograms", "He3Forward_Debug"}, "2-PhiDistribution-AllBins").Scale(30)) << fpc;
        const string cut = "set title 'Identified as 3He " + title + "'";
        SP2((HS == DATA) ? "SP2-data-cut1" : "").Distr(Hist2d(HS, R, {"Histograms", "He3Forward_Reconstruction"}, string("4-GeomCut-FTH1-vs-FRH1")).Scale(3, 3)) << cut << "set zrange [0:]";
        THD().Hist(Hist(HS, R, {"Histograms", "He3Forward_Debug"}, "4-PhiDistribution-AllBins").Scale(30)) << cut;
    };
    MakePlots(MC, "He3eta-gg", 500000);
    MakePlots(MC, "He3pi0pi0", 500000);
    MakePlots(MC, "He3pi0pi0pi0", 500000);
    MakePlots(DATA, "All", 70000);
    MakePlots(MC, "He3pi0", 500000);
}
