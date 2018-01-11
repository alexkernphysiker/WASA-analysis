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
#include <Genetic/fit.h>
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
        res += pp(Pt.Transform(pt.Beta()).P().M());
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
    RANDOM random;
    const string ppn_reaction = "ppn_qf_",pd_reaction="pd_";
    const auto runs = PresentRuns("E");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    const string th1 = "'Theta_1, deg'", th2 = "'Theta_2, deg'",
            e1 = "'Edep_1, GeV'", e2 = "'Edep_2, GeV'",
            planarity = "'Phi_1-Phi_2, deg'";
    const hist<> norm = Hist(MC, ppn_reaction, {"Histograms", "elastic"}, "0-Reference");
    const hist<> norm_pd = Hist(MC, pd_reaction, {"Histograms", "elastic"}, "0-Reference");
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "luminosity-central");
    Plot("ppn-copl-mc")
    .Hist(Hist(MC, ppn_reaction, {"Histograms", "elastic"}, "pair_phi_diff_0") / norm.TotalSum().val(), "ppn_{sp}")
            << "set key on" << "set title 'Coplanarity. MC'" << "set xrange [0:360]"
            << "set yrange [0:]" << "set xlabel " + planarity;
    Plot()
    .Hist(Hist(MC, ppn_reaction, {"Histograms", "elastic"}, "pair_phi_diff_1") / norm.TotalSum().val(), "ppn_{sp}")
            << "set key on" << "set title 'Coplanarity. MC. Cut'" << "set xrange [0:360]"
            << "set yrange [0:]" << "set xlabel " + planarity;
    Plot("ppn-copl-data")
    .Hist(Hist(DATA, "E", {"Histograms", "elastic"}, "pair_phi_diff_1"))
            << "set title 'Coplanarity. Data " + runmsg + "'" << "set xrange [0:360]"
            << "set yrange [0:]" << "set xlabel " + planarity;

    PlotHist2d(sp2).Distr(Hist2d(MC, pd_reaction, {"Histograms", "elastic"}, "t_vs_t_1"))
            << "set xrange [23:50]"<< "set yrange [23:80]"
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d(sp2).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_t_1"))
            << "set xrange [23:50]"<< "set yrange [23:80]"
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d(sp2).Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_t_1"))
            << "set xrange [23:50]"<< "set yrange [23:80]"
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + th2;


    PlotHist2d(sp2, "pd-tvt-mc-0").Distr(Hist2d(MC, pd_reaction, {"Histograms", "elastic"}, "t_vs_t_21"))
            << "set xrange [23:50]"<< "set yrange [23:80]"
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d(sp2, "ppn-tvt-mc-0").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_t_21"))
            << "set xrange [23:50]"<< "set yrange [23:80]"
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d(sp2, "ppn-tvt-data-0").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_t_21"))
            << "set xrange [23:50]"<< "set yrange [23:80]"
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d(sp2).Distr(Hist2d(MC, pd_reaction, {"Histograms", "elastic"}, "t_vs_e_21"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d(sp2).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "t_vs_e_21"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d(sp2).Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "t_vs_e_21"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + e1;
    PlotHist2d(sp2, "pd-eve-mc-0").Distr(Hist2d(MC, pd_reaction, {"Histograms", "elastic"}, "e_vs_e_21"))
            << "set zrange [0:]" << "set title 'MC pd'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d(sp2, "ppn-eve-mc-0").Distr(Hist2d(MC, ppn_reaction, {"Histograms", "elastic"}, "e_vs_e_21"))
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + e1 << "set ylabel " + e2;
    PlotHist2d(sp2, "ppn-eve-data-0").Distr(Hist2d(DATA, "E", {"Histograms", "elastic"}, "e_vs_e_21"))
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + e1 << "set ylabel " + e2;


    hist<> acceptance, acceptance_pd,events,events2,data_chi_sq,chi_sq_total;
    const auto diff_cs = ReadCrossSection();
    const auto p_cs = IntegrateCrossSection(diff_cs);
    Plot("pp-integrated").Line(p_cs) << "set title 'pp->pp'";
    const auto ppn_cs = pp2ppn(p_cs)*0.96;
    const auto SIGMA = ConvertCrossSections(ppn_cs);
    Plot("ppn-integrated").Line(ppn_cs.XRange(p_beam_low, p_beam_hi)) << "set title 'pd->pp+n_{sp}'";
    Plot("ppn-sigma").Hist(SIGMA)
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

        const hist<> data_copl=
            Hist(DATA,"E",{"Histograms","elastic"},string("pair_phi_diff_21-Bin-") + to_string(bin_num))
            .Scale(4).XRange(90,270);
        const hist<> data_copl_mc=
            Hist(MC,ppn_reaction,{"Histograms","elastic"},string("pair_phi_diff_21-Bin-") + to_string(bin_num))
            .Scale(4).XRange(90,270);
        const hist<> data_copl_mc2=
            Hist(MC,pd_reaction,{"Histograms","elastic"},string("pair_phi_diff_21-Bin-") + to_string(bin_num))
            .Scale(4).XRange(90,270);
        const hist<> data_copl_inside=data_copl.XRange(140,220);
        const hist<> data_copl_outside=data_copl.XExclude(140,220);
        cout << endl << Qmsg << endl;
        cout << endl;

        const auto epsilon = std_error(data_copl_mc.TotalSum().val())/N;
        const auto epsilon2 = std_error(data_copl_mc2.TotalSum().val())/N_pd;
        acceptance << make_point(Q, epsilon);
        acceptance_pd << make_point(Q, epsilon2);

        Fit<DifferentialMutations<Uncertainty>> fit(
            make_shared<FitPoints>()<<data_copl_outside,
            [](const ParamSet&X,const ParamSet&P){return Polynom<2>(X[0],P);}
        );
        fit.SetUncertaintyCalcDeltas({0.001,0.001,0.001});
        fit.Init(500,make_shared<InitialDistributions>()
            <<make_shared<DistribGauss>(10000,10000)
            <<make_shared<DistribGauss>(0,1)
            <<make_shared<DistribGauss>(-2,1)
            ,random
        );
        while(!fit.AbsoluteOptimalityExitCondition(0.0000001))fit.Iterate(random);
        cout << "Fitting: " << fit.iteration_count() << " iterations; "
            << fit.Optimality() << "<chi^2<"
            << fit.Optimality(fit.PopulationSize() - 1)
            << endl;
        const auto &P = fit.ParametersWithUncertainties();
        const auto BG=[&P](const value<>&x){return Polynom<2>(x,P);};
        data_chi_sq << point<value<>>(Q, fit.Optimality() / (data_copl_outside.size() - fit.ParamCount()));
        const hist<> data_copl_fg=data_copl_inside-BG,data_copl_bg=(data_copl*0.)+BG;
        const auto ev=data_copl_fg.TotalSum();
        const auto total_line=hist<>((data_copl_mc*ev/(N*epsilon))+data_copl_bg).toLine();
        Plot(Q.Contains(21) ? "ppn-above-data-copl" : (Q.Contains(-39) ? "ppn-below-data-copl" : ""))
            .Hist(data_copl_inside).Hist(data_copl_outside)
            .Line(data_copl_bg.toLine(),"BG from fit")
            .Line(total_line,"MC+BG")
                << "set title 'Coplanarity. Data " + runmsg+ "; "+Qmsg + "'" <<"set key on"
                << "set yrange [0:]" << "set xlabel " + planarity;
        double total_chisq=0;
        for(size_t i=0;i<data_copl.size();i++)total_chisq+=data_copl[i].Y().NumCompare(total_line[i].Y());
        chi_sq_total << make_point(Q, total_chisq / (data_copl.size()-fit.ParamCount()));
        events<<make_point(Q,ev);

        const hist<> data_time=
            Hist(DATA,"E",{"Histograms","elastic"},string("pair_time_diff_21-Bin-") + to_string(bin_num))
            .Scale(1).XRange(-50,50);
        const hist<> data_time_mc=
            Hist(MC,ppn_reaction,{"Histograms","elastic"},string("pair_time_diff_21-Bin-") + to_string(bin_num))
            .Scale(1).XRange(-50,50);
        const hist<> data_time_mc2=
            Hist(MC,pd_reaction,{"Histograms","elastic"},string("pair_time_diff_21-Bin-") + to_string(bin_num))
            .Scale(1).XRange(-50,50);
        const hist<> time_fg=data_time.XRange(-20,-7);
        const hist<> time_bg=time_fg.CloneEmptyBins().Transform([&time_fg](const value<>&x,const value<>&){
            const auto& a=time_fg.left().X();
            const auto& b=time_fg.right().X();
            const auto&xa=time_fg.left().Y();
            const auto&xb=time_fg.right().Y();
            return xa+(xb-xa)*(x-a)/(b-a);
        });
        const auto ev2=hist<>(time_fg-time_bg).TotalSum();

        Plot(Q.Contains(21) ? "ppn-above-data-time" : (Q.Contains(-39) ? "ppn-below-data-time" : ""))
            .Hist(data_time,"Data").Hist(time_bg,"Background")
                << "set title 'Time difference. Data " + runmsg+ "; "+Qmsg + "'" <<"set key on"
                << "set yrange [0:]" << "set xlabel 'time difference, ns'";
        events2<<make_point(Q,ev2);

    }
    Plot("ppn-acceptance")
        .Hist(acceptance, "ppn_{sp}").Hist(acceptance_pd, "pd")
            << "set key on" << "set title 'Acceptance'" << "set yrange [0:]" 
            << "set xlabel 'Q, MeV'" << "set ylabel 'Acceptance, n.d.'";
    Plot("ppn-chisq")
        .Hist(data_chi_sq,"BG").Hist(chi_sq_total,"Total")
            << "set xlabel 'Q, MeV'"<<"set key on"
            << "set ylabel 'chi^2/d, n.d.'"
            << "set yrange [0:]" << "unset log y";

    Plot("ppn-events")
        .Hist(events,"phi").Hist(events2,"time")
            << "set key on" << "set title 'True events count "+runmsg+"'" << "set yrange [0:]" 
            << "set xlabel 'Q, MeV'" << "set ylabel 'count, n.d.'";
    const hist<> luminosity=(events*double(trigger_elastic1.scaling)/acceptance/SIGMA);
    const hist<> luminosity2=(events2*double(trigger_elastic1.scaling)/acceptance/SIGMA);
    /*
    const auto sasha4d=Plotter::Instance().GetPoints<double>("luminosity_Q");
    Chain<point<>> machine4d;
    for(const auto&p:sasha4d)machine4d.push_back(make_point(-72.5+(p.X()*2.5),p.Y()));
    */
    const hist<> prev_luminosity = Plotter::Instance().GetPoints<value<>>("LUMINOSITYf");
    Plot("luminosity-compare")
        .Hist(luminosity, "ppn_{sp} (angular)")
        .Hist(luminosity2, "ppn_{sp} (time)", "LUMINOSITYc")
        .Hist(prev_luminosity, "3He+eta")
        //.Line(machine4d,"Sasha")
            << "set title 'Integrated luminosity (" + runmsg + ")'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:]";
    cout<<"Full luminosity estimation (phi): "<<luminosity.TotalSum()<<endl;
    cout<<"Full luminosity estimation (time): "<<luminosity2.TotalSum()<<endl;
}
