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
#include <Genetic/equation2.h>
#include <Genetic/initialconditions.h>
#include <Genetic/parabolic.h>
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
    BiSortedPoints<> result(ChainWithCount(91, 0., PI<>() / 2.0), ChainWithCount(13, 1.000, 2.200));
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
        result << make_point(P, IntT[IntT.size() - 1].Y() * PI() * 2.0);
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
    .Hist(Hist(DATA, "E", {"Histograms", "elastic"}, "pair_phi_diff_24-AllBins"))
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
    PlotHist2d<>(sp2, "pd-eve-mc-0").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "e_vs_e_22"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2, "ppn-eve-mc-0").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "e_vs_e_22"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2, "ppn-eve-data-0").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "e_vs_e_22"))
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
    PlotHist2d<>(sp2, "pd-eve-mc-1").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "e_vs_e_23"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2, "ppn-eve-mc-1").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "e_vs_e_23"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2, "ppn-eve-data-1").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "e_vs_e_23"))
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
    PlotHist2d<>(sp2, "pd-eve-mc-2").Distr(Hist2d(MC, "pd", {"Histograms", "elastic"}, "e_vs_e_24"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2, "ppn-eve-mc-2").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "e_vs_e_24"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d<>(sp2, "ppn-eve-data-2").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "e_vs_e_24"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + e1 << "set ylabel " + e2;

    Plot<>("ppn-theta2-mc")
    .Line(Hist(MC, "pd", {"Histograms", "elastic"}, "theta_s_24-AllBins").toLine() / norm_pd.TotalSum().val(), "pd")
    .Line(Hist(MC, ppn_reaction, {"Histograms", "elastic"}, "theta_s_24-AllBins").toLine() / norm.TotalSum().val(), "ppn_{sp}")
            << "set title 'MC'" << "set key on"
            << "set yrange [0:]" << "set xlabel " + th2
            << "set ylabel 'counts normalized'";
    Plot<>("ppn-theta2-data")
    .Hist(Hist(DATA, "E", {"Histograms", "elastic"}, "theta_s_24-AllBins"))
            << "set title 'Data " + runmsg + "'" << "set key on" << "set yrange [0:]" << "set xlabel " + th2;


    hist<> acceptance_l,acceptance_r, acceptance_pd_l,acceptance_pd_r,eventsl,eventsr,events1,evl,evr,events2,chi_sq;
    const auto diff_cs = ReadCrossSection();
    const auto p_cs = IntegrateCrossSection(diff_cs);
    Plot<>("pp-integrated").Line(p_cs) << "set title 'pp->pp'";
    const auto ppn_cs = pp2ppn(p_cs);
    const auto SIGMA = ConvertCrossSections(ppn_cs);
    Plot<>("ppn-integrated").Line(ppn_cs.XRange(p_beam_low, p_beam_hi)) << "set title 'pd->pp+n_{sp}'";
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
            static_cast<stringstream &>(
                stringstream()
                << "Q in [" << setprecision(3)
                << Q.min() << "; " << Q.max() << "] MeV"
            ).str();

        const hist<> mc_ppn =
            Hist(MC, ppn_reaction, {"Histograms", "elastic"}, string("theta_s_24-Bin-") + to_string(bin_num))
            .Scale(2).XRange(20, 110);
        const hist<> mc_pd =
            Hist(MC, "pd", {"Histograms", "elastic"}, string("theta_s_24-Bin-") + to_string(bin_num))
            .Scale(2).XRange(20, 110);
        const hist<> data =
            Hist(DATA, "E", {"Histograms", "elastic"}, string("theta_s_24-Bin-") + to_string(bin_num))
            .Scale(2).XRange(20, 110);

        const hist<> mc_ppn_l = mc_ppn.XRange(0,80);
        const hist<> mc_pd_l = mc_pd.XRange(0,80);
        const hist<> data_l = data.XRange(0,80);
        const hist<> mc_ppn_r = mc_ppn.XRange(80,180);
        const hist<> mc_pd_r = mc_pd.XRange(80,180);
        const hist<> data_r = data.XRange(80,180);

        const auto epsilon_l = std_error(mc_ppn_l.TotalSum().val())/N;
        const auto epsilon2_l = std_error(mc_pd_l.TotalSum().val())/N_pd;
        acceptance_l << make_point(Q, epsilon_l);
        acceptance_pd_l << make_point(Q, epsilon2_l);
        const auto epsilon_r = std_error(mc_ppn_r.TotalSum().val())/N;
        const auto epsilon2_r = std_error(mc_pd_r.TotalSum().val())/N_pd;
        acceptance_r << make_point(Q, epsilon_r);
        acceptance_pd_r << make_point(Q, epsilon2_r);

        const auto events_l = std_error(data_l.TotalSum().val());
        const auto events_r = std_error(data_r.TotalSum().val());

        cout << endl << Qmsg << endl;
        Plot<>(Q.Contains(21) ? "ppn-above-mc" : (Q.Contains(-39) ? "ppn-below-mc" : ""))
        .Hist(mc_ppn / N, "ppn_{sp}")
        .Line(mc_pd.toLine()/N_pd.val(),"pd left").Line(mc_pd_r.toLine()/N_pd.val(),"pd right")
                << "set key on" << "set title 'MC " + Qmsg + "'" << "set yrange [0:]"
                << "set xlabel " + th2 << "set ylabel 'counts normalized'";
        Plot<>(Q.Contains(21) ? "ppn-above-data" : (Q.Contains(-39) ? "ppn-below-data" : ""))
        .Hist(data_l).Hist(data_r)
                << "set key on" << "set title 'Data " + Qmsg+", "+runmsg + "'" << "set yrange [0:]"
                << "set xlabel " + th2 << "set ylabel 'counts normalized'";
        cout << endl;

        RANDOM r_eng;
        InexactEquationSolver<DifferentialMutations<Uncertainty>> solver{
            {.left=[&epsilon_l,&epsilon2_l](const ParamSet&E){return epsilon_l*E[0]+epsilon2_l*E[1];},.right=events_l},
            {.left=[&epsilon_r,&epsilon2_r](const ParamSet&E){return epsilon_r*E[0]+epsilon2_r*E[1];},.right=events_r}
        };
        const auto expected_count=(events_l+events_r).val();
        solver.SetFilter([](const ParamSet&E){return (E[0]>0)&&(E[1]>0);})
        .Init(
            200,
            make_shared<InitialDistributions>()
                << make_shared<DistribUniform>(0, expected_count*50.)
                << make_shared<DistribUniform>(0, expected_count*5),
            r_eng
        );
        while (!solver.AbsoluteOptimalityExitCondition(0.0000000001))
            solver.Iterate(r_eng);
        solver.SetUncertaintyCalcDeltas({0.01, 0.01});
        const auto &E = solver.ParametersWithUncertainties();
        chi_sq<<make_point(Q,solver.Optimality());

        evl<<make_point(Q,solver.equations()[0].left(solver.Parameters()));
        evr<<make_point(Q,solver.equations()[1].left(solver.Parameters()));
        eventsl<<make_point(Q,solver.equations()[0].right);
        eventsr<<make_point(Q,solver.equations()[1].right);

        events1<<make_point(Q,E[0]* double(trigger_elastic1.scaling));
        events2<<make_point(Q,E[1]* double(trigger_elastic1.scaling));
    }
    Plot<>("ppn-acceptance")
    .Line(acceptance_l.toLine(), "ppn_{sp} left").Line(acceptance_pd_l.toLine(), "pd left")
    .Line(acceptance_r.toLine(), "ppn_{sp} right").Line(acceptance_pd_r.toLine(), "pd right")
            << "set key on" << "set title 'Acceptance'" << "set yrange [0:]" << "set xlabel 'Q, MeV'" << "set ylabel 'Acceptance, n.d.'";
    Plot<>("ppn-chisq")
    .Hist(chi_sq)
            << "set key on" << "set title 'chi^2'" << "set yrange [0:]" << "set xlabel 'Q, MeV'" << "set ylabel 'chi^2, n.d.'";
    Plot<>("ppn-events-lr")
    .Hist(eventsl, "left").Hist(eventsr, "right")
    .Line(evl.toLine()).Line(evr.toLine())
            << "set key on" << "set title 'Registered events count'" << "set yrange [0:]" << "set xlabel 'Q, MeV'" << "set ylabel 'count, n.d.'";
    Plot<>("ppn-events")
    .Hist(events1, "ppn_{sp}").Hist(events2, "pd")
            << "set key on" << "set title 'True events count'" << "set yrange [0:]" << "set xlabel 'Q, MeV'" << "set ylabel 'count, n.d.'";
    const auto luminosity=events1/SIGMA;
    Plot<>("ppn-luminosity").Hist(luminosity)
            << "set title 'Integrated luminosity (" + runmsg + ")'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:]";
    const hist<> prev_luminosity = Plotter<>::Instance().GetPoints<value<>>("LUMINOSITYf");
    const hist<> estimate_full_luminosity = luminosity * runs.second / runs.first;
    Plot<>("luminosity-compare")
    .Hist(estimate_full_luminosity, "ppn_{sp}", "LUMINOSITYc")
    .Hist(prev_luminosity, "3He+eta")
            << "set title 'Integrated luminosity estimation for all runs'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:]";
}


