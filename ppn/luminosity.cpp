// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <fstream>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <math_h/vectors.h>
#include <Genetic/searchmin.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Genetic/parabolic.h>
#include <Genetic/fit.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Kinematics/reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
const BiSortedPoints<> ReadCrossSection()
{
    BiSortedPoints<> result(ChainWithCount(91, 0., PI<>()/2.0), ChainWithCount(13, 1.000, 2.200));
    for (size_t degree = 0; degree <= 90; degree++) {
        ifstream file("crosssections/Theta_" + to_string(degree) + ".txt");
        for (double E = 0, C = 0; (file >> E >> C); result.Bin(degree, (size_t(E) - 1000) / 100) = C);
        file.close();
    }
    return result;
}
const SortedPoints<> IntegrateCrossSection(const BiSortedPoints<> &angular)
{
    SortedPoints<> result;
    result << make_point(0.0, 0.0) << make_point(angular.Y()[0] - 0.001, 0.0);
    result << make_point(5.0, 0.0) << make_point(angular.Y()[angular.Y().size() - 1] + 0.001, 0.0);
    for (size_t index = 0; index < angular.Y().size(); index++) {
        const auto &P = angular.Y()[index];
        const auto IntT = Int_Trapez_Table(angular.CutX(index)
        * [](const double & th) {
            return sin(th);
        });
        result << make_point(P, IntT[IntT.size() - 1].Y() * PI()*2.0);
    }
    return result;
}
const Points<> ReadPf()
{
    Points<> result;
    ifstream file("crosssections/pfermi.txt");
    for (double P = 0, D = 0; (file >> P >> D); result.push_back(make_point(P, D)));
    file.close();
    return result;
}
const double Calculate_pp2ppn(const double &pbeam, const function<double(double)> &pp)
{
    static const RandomValueTableDistr<> PF = ReadPf();
    RANDOM R;
    double res = 0;
    const size_t count = 10000;
    const auto Pt = lorentz_byPM(Z<>() * pbeam, Particle::p().mass()),
               T = lorentz_byPM(Zero<>(), Particle::d().mass());
    for (size_t i = 0; i < count; i++) {
        const auto
        nt = lorentz_byPM(randomIsotropic<3>(R) * PF(R), Particle::n().mass()),
        pt = T - nt;
        res += pp(Pt.Transform(pt.Beta()).S().M());
    }
    return res / count;
}
const SortedPoints<> pp2ppn(const LinearInterpolation<> &pp_cs)
{
    SortedPoints<> result;
    for (double p = 1.0; p <= 2.1; p += 0.001) {
        result << make_point(p, Calculate_pp2ppn(p, pp_cs.func()));
    }
    return result;
}
const SortedPoints<value<>> ConvertCrossSections(const SortedPoints<> &momentum)
{
    LinearInterpolation<value<>> q_form;
    static const Reaction he3eta(Particle::p(), Particle::d(), {Particle::he3(), Particle::eta()});
    for (const auto &P : momentum) {
        q_form << point<value<>>(he3eta.P2Q(P.X()) * 1000., P.Y() * 1000000);
    }
    return SortedPoints<value<>>(q_form.func(), BinsByStep(-70., 2.5, +30.));
}
int main()
{
    const string ppn_reaction = "ppn_qf_";
    const auto runs = PresentRuns("E");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    const string th1 = "'Theta_1, deg'", th2 = "'Theta_2, deg'", e1 = "'Edep_1, GeV'", e2 = "'Edep_2, GeV'",
                 thth = "'Theta_1+1.6Theta_2, deg'", planarity = "'|Phi_1-Phi_2-180^o|, deg'";
    const hist<> norm = Hist(MC, ppn_reaction, {"Histograms", "elastic"}, "0-Reference");
    const hist<> norm_pd = Hist(MC, "pd", {"Histograms", "elastic"}, "0-Reference");
    Plotter<>::Instance().SetOutput(ENV(OUTPUT_PLOTS), "luminosity-central");
    Plot<>("ppn-copl-mc")
    .Hist(Hist(MC, "pd", {"Histograms", "elastic"}, "pair_phi_diff_0") / norm_pd.TotalSum().val(), "pd")
    .Hist(Hist(MC, ppn_reaction, {"Histograms", "elastic"}, "pair_phi_diff_0") / norm.TotalSum().val(), "ppn_{sp}")
            << "set key on" << "set title 'Coplanarity. MC'" << "set yrange [0:]" << "set xlabel " + planarity;
    Plot<>()
    .Hist(Hist(MC, "pd", {"Histograms", "elastic"}, "pair_phi_diff_1") / norm_pd.TotalSum().val(), "pd")
    .Hist(Hist(MC, ppn_reaction, {"Histograms", "elastic"}, "pair_phi_diff_1") / norm.TotalSum().val(), "ppn_{sp}")
            << "set key on" << "set title 'Coplanarity. MC. Cut'" << "set yrange [0:]" << "set xlabel " + planarity;
    Plot<>("ppn-copl-data")
    .Hist(Hist(DATA, "E", {"Histograms", "elastic"}, "pair_phi_diff_0"))
    .Hist(Hist(DATA, "E", {"Histograms", "elastic"}, "pair_phi_diff_1"))
    .Hist(Hist(DATA, "E", {"Histograms", "elastic"}, "pair_phi_diff_24"))
            << "set title 'Coplanarity. Data " + runmsg + "'" << "set yrange [0:]" << "set xlabel " + planarity;

    PlotHist2d<>(sp2).Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "t_vs_t_1"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_t_1"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_t_1"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + th2;


    PlotHist2d<>(sp2, "pd-tvt-mc-0").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "t_vs_t_22"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2, "ppn-tvt-mc-0").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_t_22"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2, "ppn-tvt-data-0").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_t_22"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "t_vs_e_22"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_e_22"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_e_22"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2,"pd-eve-mc-0").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "e_vs_e_22"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2,"ppn-eve-mc-0").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "e_vs_e_22"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2,"ppn-eve-data-0").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "e_vs_e_22"))    
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + e1 << "set ylabel " + e2;
            
    PlotHist2d<>(sp2, "pd-tvt-mc-1").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "t_vs_t_23"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2, "ppn-tvt-mc-1").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_t_23"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2, "ppn-tvt-data-1").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_t_23"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "t_vs_e_23"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_e_23"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_e_23"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2,"pd-eve-mc-1").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "e_vs_e_23"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2,"ppn-eve-mc-1").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "e_vs_e_23"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2,"ppn-eve-data-1").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "e_vs_e_23"))    
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + e1 << "set ylabel " + e2;
            
    PlotHist2d<>(sp2, "pd-tvt-mc-2").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "t_vs_t_24"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2, "ppn-tvt-mc-2").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_t_24"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2, "ppn-tvt-data-2").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_t_24"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "t_vs_e_24"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_e_24"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_e_24"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d<>(sp2,"pd-eve-mc-2").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "e_vs_e_24"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2,"ppn-eve-mc-2").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "e_vs_e_24"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2,"ppn-eve-data-2").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "e_vs_e_24"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + e1 << "set ylabel " + e2;
            
    Plot<>("ppn-sumofthetas-mc")
    .Line(Hist(MC, "pd", {"Histograms", "elastic"}, "theta_sum_24-AllBins").toLine() / norm_pd.TotalSum().val(), "pd")
    .Line(Hist(MC, ppn_reaction, {"Histograms", "elastic"}, "theta_sum_24-AllBins").toLine() / norm.TotalSum().val(), "ppn_{sp}")
            << "set title 'MC'" << "set key on"
            << "set yrange [0:]" << "set xlabel " + thth
            << "set ylabel 'counts normalized'";
    Plot<>("ppn-sumofthetas-data")
    .Hist(Hist(DATA, "E", {"Histograms", "elastic"}, "theta_sum_24-AllBins"))
            << "set title 'Data " + runmsg + "'" << "set key on" << "set yrange [0:]" << "set xlabel " + thth;


    RANDOM r_eng;
    hist<> acceptance, acceptance_pd, chi_sq, luminosity, el_cs;
    vector<hist<>> fit_params;
    const auto diff_cs = ReadCrossSection();
    const auto p_cs = IntegrateCrossSection(diff_cs);
    Plot<>("pp-integrated").Line(p_cs) << "set title 'pp->pp'";
    const auto ppn_cs = pp2ppn(p_cs);
    const auto SIGMA = ConvertCrossSections(ppn_cs);
    Plot<>("ppn-integrated").Line(ppn_cs.XRange(p_beam_low,p_beam_hi)) << "set title 'pd->pp+n_{sp}'";
    Plot<>("ppn-sigma").Hist(SIGMA)
            << "set title 'ppn_{sp} cross section'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'cross section, nb'"
            << "set xrange [-70:30]" << "set yrange [0:]";
    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const auto &N = norm[bin_num].Y();
        const auto &N_pd = norm_pd[bin_num].Y();
        const string Qmsg =
            static_cast<stringstream &>(stringstream()
                                        << "Q in [" << setprecision(3)
                                        << Q.min() << "; " << Q.max() << "] MeV"
                                       ).str();

        const hist<> mc_ppn =
            Hist(MC, ppn_reaction, {"Histograms", "elastic"}, string("theta_s_24-Bin-") + to_string(bin_num))
            .Scale(2).XRange(20, 200);
        const hist<> mc_pd =
            Hist(MC, "pd", {"Histograms", "elastic"}, string("theta_s_24-Bin-") + to_string(bin_num))
            .Scale(2).XRange(20, 200);
        const hist<> data =
            Hist(DATA, "E", {"Histograms", "elastic"}, string("theta_s_24-Bin-") + to_string(bin_num))
            .Scale(2).XRange(20, 200);
        const hist<> nmc_ppn = mc_ppn / N;
        const hist<> nmc_pd = mc_pd / N_pd;
        const auto epsilon = nmc_ppn.TotalSum();
        acceptance << make_point(Q, epsilon);
        acceptance_pd << make_point(Q, nmc_pd.TotalSum());
        cout << endl << Qmsg << endl;
        Plot<>(Q.Contains(21) ? "ppn-above-mc" : (Q.Contains(-39) ? "ppn-below-mc" : ""))
        .Hist(nmc_ppn, "ppn_{sp}").Hist(nmc_pd, "pd")
                << "set key on" << "set title 'MC " + Qmsg + "'" << "set yrange [0:]"
                << "set xlabel " + thth << "set ylabel 'counts normalized'";
        cout << endl;
        const std::function<const double(const double &, const ParamSet &)>
        BG = [](const double & x, const ParamSet & P) {
            return FermiFunc(x, P[2], P[3]) * P[4];
        };
        /*
        const auto &data_count = data.TotalSum().val();
        SearchMin<DifferentialMutations<Uncertainty>>
        FitData([&data, &nmc_ppn, &nmc_pd, BG](const ParamSet & P) {
            double res = 0;
            for (size_t i = 0; i < data.size(); i++) {
                const double x = data[i].X().val();
                const auto practic = data[i].Y();
                const auto theor =
                    nmc_ppn[i].Y() * P[0]
                    + nmc_pd[i].Y() * P[1]
                    + BG(x, P);
                res += practic.NumCompare(theor);
            }
            return res;
        });
        FitData.SetFilter([](const ParamSet & P) {
            return (P[0] > 0) && (P[1] > 0)
                   && (P[2] > 50) && (P[2] < 80)
                   && (P[3] < 0) && (P[4] > 0);
        });
        FitData.Init(200, make_shared<InitialDistributions>()
                     << make_shared<DistribUniform>(0, 10.*data_count)
                     << make_shared<DistribUniform>(0, 10.*data_count)
                     << make_shared<DistribUniform>(60, 70)
                     << make_shared<DistribUniform>(-5, 0)
                     << make_shared<DistribUniform>(0, 0.01 * data_count)
                     , r_eng
                    );
        FitData.SetUncertaintyCalcDeltas(parEq(FitData.ParamCount(), 0.1));
        while (!FitData.AbsoluteOptimalityExitCondition(0.0000000001)) FitData.Iterate(r_eng);
        cout << "DATA: " << FitData.iteration_count() << " iterations; "
             << FitData.Optimality() << "<chi^2<"
             << FitData.Optimality(FitData.PopulationSize() - 1)
             << endl;
        const auto &P = FitData.ParametersWithUncertainties();
        const auto &p = FitData.Parameters();
        for (size_t i = 0; i < P.size(); i++) {
            if (fit_params.size() == i)fit_params.push_back(hist<>());
            fit_params[i] << point<value<>>(Q, P[i]);
        }
        const SortedPoints<>
        PPN = nmc_ppn.toLine() * p[0], PD = nmc_pd.toLine() * p[1],
        BackGround = data.toLine().Clone().Transform([BG, &p](const double & x, const double &) {
            return BG(x, p);
        });
        Plot<>(Q.Contains(21) ? "ppn-above-fit" : (Q.Contains(-39) ? "ppn-below-fit" : ""))
        .Hist(data, "DATA")
        .Line(PPN + PD + BackGround, "total fit")
        .Line(PD + BackGround, "pd+background")
        .Line(BackGround, "background")
                << "set key on" << "set title 'Data " + Qmsg + "(" + runmsg + ")'" << "set yrange [0:]"
                << "set xlabel " + thth << "set ylabel 'counts'";

        chi_sq << make_point(Q, FitData.Optimality() / (data.size() - FitData.ParamCount()));
        
        const auto pd = nmc_pd * P[1],
        background = data.Clone().Transform([BG, &P](const value<> &x, const value<> &) {
            return P[4] * FermiFunc(x.val(), P[2].val(), P[3].val());
        });
        const hist<> foreground = data - pd - background;
        Plot<>(Q.Contains(21) ? "ppn-above-subtr" : (Q.Contains(-39) ? "ppn-below-subtr" : ""))
        .Hist(foreground).Line(Points<> {{50, 0}, {200, 0}})
                << "set title 'Subtracted background " + Qmsg + "(" + runmsg + ")'";
                */
        const auto L =
           data.TotalSum() / SIGMA[bin_num].Y()/epsilon
            * double(trigger_elastic1.scaling);
        //const auto EL =
         //   P[1] / L
         //   * double(trigger_elastic1.scaling);
        luminosity << make_point(Q, L);
        //el_cs << make_point(Q, EL);
    }
    Plot<>("ppn-acceptance").Hist(acceptance, "ppn_{sp}").Hist(acceptance_pd, "pd") << "set key on"
            << "set title 'Acceptance'" << "set yrange [0:]" << "set xlabel 'Q, MeV'" << "set ylabel 'Acceptance, n.d.'";
/*
    for (size_t i = 0; i < fit_params.size(); i++) {
        Plot<>().Hist(fit_params[i])
                << "set xlabel 'Q, MeV'"
                << "set ylabel 'parameter" + to_string(i) + "'";
    }

    Plot<>("ppn-chisq").Hist(chi_sq, "DATA")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'chi^2/d, n.d.'"
            << "set yrange [0:]" << "unset log y";
*/
    Plot<>("ppn-luminosity").Hist(luminosity)
            << "set title 'Integrated luminosity (" + runmsg + ")'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:]";
    const hist<> prev_luminosity = Plotter<>::Instance().GetPoints<value<>>("LUMINOSITYf");
    const hist<> estimate_full_luminosity = luminosity * runs.second / runs.first;
    Plot<>("luminosity-compare")
    .Hist(estimate_full_luminosity, "ppn_{sp}")
    .Hist(prev_luminosity, "3He+eta")
            << "set title 'Integrated luminosity estimation for all runs'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:]";
    SearchMin<DifferentialMutations<Uncertainty>> Shadowing([&estimate_full_luminosity, &prev_luminosity](const ParamSet & P) {
        double res = 0;
        for (size_t i = 0; i < prev_luminosity.size(); i++) {
            const size_t ii = estimate_full_luminosity.size() - (prev_luminosity.size() - i);
            res += (estimate_full_luminosity[ii].Y() * P[0]).NumCompare(prev_luminosity[i].Y());
        }
        return res;
    });
    Shadowing.SetFilter([](const ParamSet & P) {
        return P[0] > 0;
    });
    Shadowing.Init(10, make_shared<InitialDistributions>() << make_shared<DistribUniform>(1.0, 2.0), r_eng);
    while (!Shadowing.AbsoluteOptimalityExitCondition(0.0000000001))Shadowing.Iterate(r_eng);
    Shadowing.SetUncertaintyCalcDeltas({0.001});
    cout << "Shadowing effect coefficient: " << Shadowing.ParametersWithUncertainties()[0] << endl;
    cout << "chi^2/d: " << Shadowing.Optimality() / (prev_luminosity.size() - Shadowing.ParamCount()) << endl;
    const hist<> full_luminosity = estimate_full_luminosity * Shadowing.ParametersWithUncertainties()[0];
    Plot<>("luminosity-compare-with-shadowing")
    .Hist(full_luminosity, "ppn_{sp}", "LUMINOSITYc")
    .Hist(prev_luminosity, "3He+eta")
            << "set title 'Taking shadowing effect into account'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:]";
}


