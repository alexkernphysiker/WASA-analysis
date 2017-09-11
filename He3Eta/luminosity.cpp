// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Genetic/fit.h>
#include <Genetic/paramfunc.h>
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
    Plotter<>::Instance().SetOutput(ENV(OUTPUT_PLOTS), "luminosity-forward");
    vector<string> histpath_forward_reconstr = {"Histograms", "He3Forward_Reconstruction"};
    const auto runs = PresentRuns("F");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    const hist<> norm = Hist(MC, "He3eta", histpath_forward_reconstr, "0-Reference");
    hist<> luminosity, data_chi_sq, acceptance;
    vector<hist<>> parhists;
    RANDOM r_eng;
    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++)
        if (norm[bin_num].X() > 5.0) {
            const auto &Q = norm[bin_num].X();
            const auto &N = norm[bin_num].Y();
            const string Qmsg = static_cast<stringstream &>(stringstream()
                                << "Q in [" << setprecision(3)
                                << Q.min() << "; " << Q.max() << "] MeV"
                                                           ).str();
            const hist<> data = Hist(DATA, "F", histpath_forward_reconstr,
                                     string("MissingMass-Bin-") + to_string(bin_num)
                                    ).XRange(0.53, 0.57);
            const auto chain = ChainWithStep(0.53, 0.0001, 0.57);
            const auto cut = make_pair(0.541, 0.554);
            const hist<> mc_unnorm = Hist(MC, "He3eta", histpath_forward_reconstr,
                                          string("MissingMass-Bin-") + to_string(bin_num)
                                         ).XRange(0.53, 0.57);
            const hist<> mc = mc_unnorm / N;
            acceptance << point<value<>>(Q, mc.TotalSum());
            Plot<>(Q.Contains(21) ? "He3eta-mc" : "")
            .Hist(mc)
                    << "set key on" << "set title '" + Qmsg + ",3He+eta MC'"
                    << "set xlabel 'Missing mass, GeV'"
                    << "set ylabel 'acceptance density, GeV^{-1}'"
                    << "set yrange [0:]" << "unset log y";
            const auto &data_count = data.TotalSum().val();
            const auto BG = [&data_count](const ParamSet & X, const ParamSet & P) {
                const double res = data_count * Polynom(X[0], P, 3, 0);
                return (res > 0) ? res : 0.0;
            };
            const auto peak_reg = data.XRange(cut.first, cut.second);
            const auto data_bg = data.XExclude(cut.first, cut.second);
            cout << endl << Qmsg << endl << endl;
            Fit<AbsoluteMutations<DifferentialMutations<Uncertainty>>> FIT(make_shared<FitPoints>(data_bg), BG);
            FIT
            .SetAbsoluteMutationCoefficients({1.0, 1.0, 1.0, 1.0})
            .SetAbsoluteMutationsProbability(0.2)
            .SetUncertaintyCalcDeltas({0.1, 0.1, 0.1, 0.1})
            .SetFilter([BG, &peak_reg](const ParamSet & P) {
                bool valid = (BG({peak_reg.left().X().val()}, P) > 0.0);
                if (valid)for (const auto &p : peak_reg) {
                        valid &= (p.Y().max() > BG({p.X().val()}, P));
                    }
                return valid;
            });
            auto init = make_shared<InitialDistributions>()
                        << make_shared<DistribGauss>(0, 5)
                        << make_shared<DistribGauss>(0, 5)
                        << make_shared<DistribGauss>(0, 5)
                        << make_shared<DistribGauss>(0, 5)
                        ;
            FIT.Init(300, init, r_eng);
            cout << endl;
            while (!FIT.AbsoluteOptimalityExitCondition(0.0000001))FIT.Iterate(r_eng);
            cout << "Fitting: " << FIT.iteration_count() << " iterations; "
                 << FIT.Optimality() << "<chi^2<"
                 << FIT.Optimality(FIT.PopulationSize() - 1)
                 << endl;
            const auto &P = FIT.ParametersWithUncertainties();
            if (parhists.size() == 0) {
                for (size_t i = 0; i < P.size(); i++)
                    parhists.push_back(hist<>());
            }
            for (size_t i = 0; i < P.size(); i++)
                parhists[i] << point<value<>>(Q, P[i]);
            data_chi_sq << point<value<>>(Q, FIT.Optimality() / (data.size() - FIT.ParamCount()));
            cout << endl;
            Plot<> exp_plot(Q.Contains(21) ? "He3eta-fit" : "");
            exp_plot.Hist(data).Hist(data_bg)
                    << "set key on" << "set title '" + Qmsg + ", " + runmsg + "'"
                    << "set xlabel 'Missing mass, GeV'"
                    << "set ylabel 'counts'"
                    << "set yrange [-200:]" << "unset log y";
            const SortedPoints<>
            background([&FIT, BG](double x) {
                return BG({x}, FIT.Parameters());
            }, chain);
            exp_plot.Line(background, "background");
            hist<> bg;
            for (const auto &po : data) {
                const auto &x = po.X().val();
                const ParamSet &p = FIT.Parameters();
                double v = BG({x}, p), u = 0;
                for (size_t index = 0; index < p.size(); index++) {
                    ParamSet p1 = p, p2 = p1;
                    const auto &delta = P[index].uncertainty();
                    p1(index) += delta;
                    p2(index) -= delta;
                    u += pow(BG({x}, p1) - BG({x}, p2), 2);
                }
                bg << point<value<>>(po.X(), {v, sqrt(u)});
            }
            hist<> clean = data - bg;
            Plot<> subplot(Q.Contains(21) ? "He3eta-subtract" : "");
            subplot.Hist(clean);
            subplot.Hist(clean = clean.XRange(cut.first, cut.second)).Object("0 title \"\"")
                    << "set key on" << "set title '" + Qmsg + ", " + runmsg + "'"
                    << "set xlabel 'Missing mass, GeV'"
                    << "set ylabel 'counts'"
                    << "set yrange [-200:]" << "unset log y";

            luminosity << point<value<>>(Q,
                                         ((clean.TotalSum() / mc.TotalSum())
                                          *trigger_he3_forward.scaling
                                          / he3eta_sigma()(Q))
                                        );
        }
    for (size_t i = 0; i < parhists.size(); i++)
        Plot<>().Hist(parhists[i])
                << "set xlabel 'Q, MeV'"
                << "set ylabel 'parameter" + to_string(i) + "'";
    Plot<>("He3eta-chisq").Hist(data_chi_sq)
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'chi^2/d, n.d.'"
            << "set yrange [0:]" << "unset log y";

    Plot<>("He3eta-cross-section")
    .Hist(he3eta_sigma(), "Data from other experiments")
    .Hist(hist<>(he3eta_sigma().func(), BinsByStep(-70.0, 2.5, 30.0)), "Interpolation", "CS-He3eta-assumed")
            << "set title 'Cross section of He3eta used in the calculations'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'sigma(^3He eta), nb'"
            << "set xrange [-20:45]" << "set yrange [0:600]";

    Plot<>("He3eta-acceptance").Hist(acceptance)
            << "set title '3He+eta acceptance'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'acceptance, n.d.'"
            << "set xrange [0:30]" << "set yrange [0:1]";

    Plot<>("He3eta-luminosity").Hist(luminosity)
            << "set title 'Integrated luminosity (" + runmsg + ")'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [0:30]" << "set yrange [0:]";

    Plot<>("luminosity-compare")
    .Hist(luminosity * runs.second / runs.first, "3He+eta", "LUMINOSITYf")
    .Hist(Plotter<>::Instance().GetPoints4("LUMINOSITYc"), "ppn_{sp}")
            << "set title 'Integrated luminosity estimation for all runs'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:]";
}
