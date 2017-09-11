// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Genetic/searchmin.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Genetic/parabolic.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include "he3eta.h"
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    Plotter<>::Instance().SetOutput(ENV(OUTPUT_PLOTS), "background-reactions-forward");
    const auto runs = PresentRuns("F");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    vector<string> histpath_forward_reconstr = {"Histograms", "He3Forward_Reconstruction"};
    vector<string> reaction = {"He3eta", "He3pi0pi0pi0", "He3pi0pi0", "He3pi0"};
    vector<hist<>> norm;
    for (const string &r : reaction)
        norm.push_back(Hist(MC, r, histpath_forward_reconstr, "0-Reference"));
    vector<hist<>> acceptance, acceptance2, true_events;
    for (size_t i = 0; i < norm.size(); i++) {
        acceptance.push_back(hist<>());
        acceptance2.push_back(hist<>());
        true_events.push_back(hist<>());
    }
    hist<> events_count;
    RANDOM r_eng;
    for (size_t bin_num = 0, bin_count = norm[0].size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[0][bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
                            << "Q in [" << setprecision(3)
                            << Q.min() << "; " << Q.max() << "] MeV"
                                                       ).str();
        auto transform = [](const hist<> &h) {
            return h.XRange(0.0, 0.6);
        };
        auto transform2 = [](const hist<> &h) {
            return h.XRange(0.4, 0.6);
        };
        const hist<> data = transform(Hist(DATA, "F", histpath_forward_reconstr, string("MissingMass-Bin-") + to_string(bin_num)));
        Plot<>(Q.Contains(21) ? "He3forward-above-mm" : (Q.Contains(-39) ? "He3forward-below-mm" : "")).Hist(data)
                << "set key on" << "set title '" + Qmsg + ", " + runmsg + "'"
                << "set xlabel 'Missing mass, GeV'"
                << "set ylabel 'counts'"
                << "set yrange [0:]" << "unset log y";
        const hist<> data2 = transform2(data);
        Plot<>(Q.Contains(21) ? "He3forward-above-mm2" : (Q.Contains(-39) ? "He3forward-below-mm" : "")).Hist(data2)
                << "set key on" << "set title '" + Qmsg + ", " + runmsg + "'"
                << "set xlabel 'Missing mass, GeV'"
                << "set ylabel 'counts'"
                << "set yrange [0:]" << "unset log y";
        {
            Plot<> th_plot(Q.Contains(21) ? "He3forward-above-mc" : (Q.Contains(-39) ? "He3forward-below-mc" : "")),
                 th_plot2(Q.Contains(21) ? "He3forward-above-mc2" : (Q.Contains(-39) ? "He3forward-below-mc2" : ""));
            th_plot << "set key on" << "set xlabel 'Missing mass, MeV'"
                    << "set title '" + Qmsg + "'"
                    << "set ylabel 'acceptance density, GeV^{-1}'"
                    << "set yrange [0:20]";
            th_plot2 << "set key on" << "set xlabel 'Missing mass, MeV'"
                     << "set title '" + Qmsg + "'"
                     << "set ylabel 'acceptance density, GeV^{-1}'"
                     << "set yrange [0:20]";
            for (size_t i = 0; i < reaction.size(); i++) {
                const hist<> react_sim = transform(Hist(MC, reaction[i], histpath_forward_reconstr, string("MissingMass-Bin-") + to_string(bin_num)));
                const hist<> react_sim2 = transform2(react_sim);
                auto N = norm[i][bin_num].Y();
                th_plot.Line(react_sim.toLine() / (N.val()*react_sim[0].X().uncertainty() * 2.), reaction[i]);
                th_plot2.Line(react_sim2.toLine() / (N.val()*react_sim2[0].X().uncertainty() * 2.), reaction[i]);
                const auto MN = value<>::std_error(react_sim.TotalSum().val());
                const auto MN2 = value<>::std_error(react_sim2.TotalSum().val());
                if (N.Above(0)) {
                    acceptance[i] << point<value<>>(Q, MN / N);
                    acceptance2[i] << point<value<>>(Q, MN2 / N);
                } else {
                    acceptance[i] << point<value<>>(Q, 0.0);
                    acceptance2[i] << point<value<>>(Q, 0.0);
                }
            }
        }
        cout << endl << Qmsg << endl << endl;
        events_count << point<value<>>(Q, value<>::std_error(data2.TotalSum().val()));
    }
    Plot<> acc("He3forward-acceptance"), acc2("He3forward-acceptance2");
    acc << "set key on"
        << "set yrange [0:1.0]" << "unset log y"
        << "set xlabel 'Q, MeV'"
        << "set ylabel 'Acceptance, n.d.'";
    acc2 << "set key on"
         << "set yrange [0:1.0]" << "unset log y"
         << "set xlabel 'Q, MeV'"
         << "set ylabel 'Acceptance, n.d.'";
    for (size_t i = 0; i < reaction.size(); i++) {
        acc.Hist(acceptance[i], reaction[i]);
        acc2.Hist(acceptance2[i], reaction[i]);
    }
    const auto luminosity = Plotter<>::Instance().GetPoints4("LUMINOSITYc");
    const auto he3etacs = Plotter<>::Instance().GetPoints4("CS-He3eta-assumed");
    const auto he3pi0pi0pi0cs = Plotter<>::Instance().GetPoints4("CS-He3pi0pi0pi0");
    const hist<> he3eta_events = luminosity * (runs.first / runs.second) / double(trigger_he3_forward.scaling) * (he3etacs * acceptance2[0]);
    const hist<> he3pi0pi0pi0_events = luminosity * (runs.first / runs.second) / double(trigger_he3_forward.scaling) * (he3pi0pi0pi0cs * acceptance2[1]);
    Plot<>("He3forward-events")
    .Hist(he3eta_events, "3He+eta")
    .Hist(he3eta_events + he3pi0pi0pi0_events, "3He+eta and 3He+3pi0")
    .Hist(events_count, "data")
            << "set key on"
            << "set yrange [0:]" << "unset log y"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'events count'";
}
