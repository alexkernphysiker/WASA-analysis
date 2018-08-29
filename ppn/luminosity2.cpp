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
#include <math_h/sigma3.h>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
#include <Genetic/uncertainties.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Parameters/parameters.h>
#include <Parameters/systematic.h>
#include <Kinematics/reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
const BiSortedPoints<> ReadCrossSection()
{
    BiSortedPoints<> result(ChainWithCount(91, 0., PI<>() / 2.0), ChainWithStep(0.700,0.050,3.000));
    for (size_t degree = 0; degree <= 90; degree++) {
        ifstream file("crosssections/Theta_" + to_string(degree) + ".txt");
        for (double E = 0, C = 0; (file >> E >> C); result.Bin(degree, (size_t(E) - 700) / 50) = C);
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
    double res = 0;
    const size_t count = 10000;
    const auto Pt = lorentz_byPM(Z<>() * pbeam, Particle::p().mass()),
               T = lorentz_byPM(Zero<>(), Particle::d().mass());
    for (size_t i = 0; i < count; i++) {
        const auto
        nt = lorentz_byPM(randomIsotropic<3>() * PF(), Particle::n().mass()),
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
    const string ppn_reaction = "ppn_qf_";
    const auto runs = PresentRuns("All");
    const string runmsg = (runs.first==runs.second)?"":"("+to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs)";
    const string th1 = "'Theta_1, deg'", th2 = "'Theta_2, deg'",
            e1 = "'Edep_1, GeV'", e2 = "'Edep_2, GeV'",
            planarity = "'Phi_1-Phi_2, deg'";
    const hist<> norm = Hist(MC, ppn_reaction, {"Histograms", "quasielastic"}, "0-Reference");
    cout << "2D plots"<<endl;
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "luminosity-central-v2");

    PlotHist2d(sp2, "ppn-v2-trackid-mc-1",5).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "quasielastic"}, "track_id-2"))
	     << "set xlabel 'E Calorimeter, GeV'"<< "set ylabel 'E PSB, GeV'"
            << "set xrange [0:0.5]"<< "set yrange [0:0.005]"<<"set xtics 0.1"
            << "set zrange [0:]" << "set title 'MC ppn'";
    PlotHist2d(sp2, "ppn-v2-trackid-data-1",4).Distr(Hist2d(DATA, "All", {"Histograms", "quasielastic"}, "track_id-2"))
	     << "set xlabel 'E Calorimeter, GeV'"<< "set ylabel 'E PSB, GeV'"
            << "set xrange [0:0.5]"<< "set yrange [0:0.005]"<<"set xtics 0.1"
            << "set zrange [0:]" << "set title 'Data "+runmsg+"'";

    Plot("ppn-v2-copl-mc",5)
    .Hist(Hist(MC, ppn_reaction, {"Histograms", "quasielastic"}, "pair_phi_diff_0-AllBins") / norm.TotalSum().val())
            << "set key on" << "set title 'MC'" << "set xrange [90:270]"<<"set xtics 30"
            << "set yrange [0:]" << "set xlabel " + planarity;
    Plot("ppn-v2-copl-data",5)
    .Hist(Hist(DATA, "All", {"Histograms", "quasielastic"}, "pair_phi_diff_0-AllBins"))
            << "set title 'Data " + runmsg + "'" << "set xrange [90:270]"<<"set xtics 30"
            << "set yrange [0:]" << "set xlabel " + planarity;

    Plot("ppn-v2-dt-mc",5)
    .Hist(Hist(MC, ppn_reaction, {"Histograms", "quasielastic"}, "pair_time_diff_0-AllBins") / norm.TotalSum().val(), "ppn_{sp}")
            << "set key on" << "set title 'MC'" << "set xrange [-25:5]"<< "set yrange [0:]" <<"set xtics 5"
	    << "set xlabel 'Time difference, ns'"<<"set ylabel 'Events count'";
    Plot("ppn-v2-dt-data",5)
        .Hist(Hist(DATA, "All", {"Histograms", "quasielastic"}, "pair_time_diff_0-AllBins")/1000.,"All")
        .Hist(Hist(DATA, "All", {"Histograms", "quasielastic"}, "pair_time_diff_1-AllBins")/1000.,"theta cut")
        .Line(toLine(Hist(DATA, "All", {"Histograms", "quasielastic"}, "pair_time_diff_0-AllBins"))/1000.)
        .Line(toLine(Hist(DATA, "All", {"Histograms", "quasielastic"}, "pair_time_diff_1-AllBins"))/1000.)
        .Line(Points<>{
            {getParameter(ppn_t1), 1200.},{getParameter(ppn_t1), 0.0},
            {getParameter(ppn_t2), 0.0},{getParameter(ppn_t2), 1200.}
        },"time cut")<<"set xtics 5"
            << "set title '" + runmsg + "'" 
	    << "set xrange [-25:5]"<< "set yrange [0:]"<<"set key on"
	    << "set xlabel 'Time difference, ns'"
	    << "set ylabel 'Events count, 10^3'";

    Plot("ppn-v2-mm-mc",5)
        .Hist(Hist(MC, ppn_reaction, {"Histograms", "quasielastic"}, "pp_mm_3") / norm.TotalSum().val())
            << "set key on" << "set title 'Missing mass'" << "set yrange [0:]"<<"set xrange [0.5:1.5]"<<"set xtics 0.1" ;
    Plot("ppn-v2-mm-data",5)
        .Hist(Hist(DATA, "All", {"Histograms", "quasielastic"}, "pp_mm_3"))
            << "set key on" << "set title 'Missing mass'" << "set yrange [0:]"<<"set xrange [0.5:1.5]"<<"set xtics 0.1" ;


    PlotHist2d(sp2, "ppn-v2-tvt-mc-1",5).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "quasielastic"}, "t_vs_t_1"))
            << "set xrange [23:40]"<< "set yrange [30:70]"<<"set xtics 5"
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d(sp2, "ppn-v2-tvt-data-1",5).Distr(Hist2d(DATA, "All", {"Histograms", "quasielastic"}, "t_vs_t_1"))
            << "set xrange [23:40]"<< "set yrange [30:70]"<<"set xtics 5"
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d(sp2, "ppn-v2-tvt-mc-2",5).Distr(Hist2d(MC, ppn_reaction, {"Histograms", "quasielastic"}, "t_vs_t_3"))
            << "set xrange [23:40]"<< "set yrange [30:70]"<<"set xtics 5"
            << "set zrange [0:]" << "set title 'MC ppn_{sp}'" << "set xlabel " + th1 << "set ylabel " + th2;
    PlotHist2d(sp2, "ppn-v2-tvt-data-2",5).Distr(Hist2d(DATA, "All", {"Histograms", "quasielastic"}, "t_vs_t_3"))
            << "set xrange [23:40]"<< "set yrange [30:70]"<<"set xtics 5"
            << "set zrange [0:]" << "set title 'Data " + runmsg + "'" << "set xlabel " + th1 << "set ylabel " + th2;

    ext_hist<2> lum_bc_z,lum_bc_p,lum_bc_m,acceptance,luminosity;
    cout << "Cross sections"<<endl;
    const auto diff_cs = ReadCrossSection();
    const auto p_cs = IntegrateCrossSection(diff_cs);
    Plot("pp-v2-integrated",4).Line(p_cs) << "set title 'pp->pp'"<<"set xtics 10";
    const auto ppn_cs = pp2ppn(p_cs)*0.955;//shading effect
    const auto SIGMA = ConvertCrossSections(ppn_cs);
    Plot("ppn-v2-integrated",4).Line(ppn_cs.XRange(p_beam_low, p_beam_hi)) << "set title 'pd->pp+n_{sp}'"<<"set xtics 10";
    Plot("ppn-v2-sigma",4).Hist(SIGMA)
        << "set title 'ppn_{sp} cross section'"
        << "set key on" << "set xlabel 'Q, MeV'"
        << "set ylabel 'cross section, nb'"<<"set xtics 10"
        << "set xrange [-70:30]" << "set yrange [0:]";
    cout << "Binning"<<endl;
    {
        for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
            const auto &Q = norm[bin_num].X();
            const auto &N = norm[bin_num].Y();
            const string Qmsg =
                static_cast<stringstream &>(
                    stringstream()
                    << "Q in [" << setprecision(3)
                    << Q.min() << "; " << Q.max() << "] MeV"
                ).str();
            size_t counter=0;
            acceptance << make_point(Q,RawSystematicError({pbeam_corr,ppn_th1,ppn_th2,ppn_t1,ppn_t2},[&N,&ppn_reaction,&bin_num](const string&suffix){//only for plot
                const hist<> data_copl_mc=
                    Hist(MC,ppn_reaction,{"Histograms","quasielastic"},string("pair_phi_diff_")+to_string(3)+string("-Bin-") + to_string(bin_num),suffix)
                    .Scale(4).XRange(90,270);
                const auto acc=data_copl_mc.TotalSum()/N;
                return extend_value<2,2>(acc);
            })());
            const auto LC=[
                    &counter,&Q,&N,&Qmsg,&bin_num,
                    &runmsg,&ppn_reaction,&planarity,
                    &SIGMA
                ](const string&suffix){
                counter++;
                cout << endl << Qmsg<<suffix << endl;
                const hist<> data_copl=
                    Hist(DATA,"All",{"Histograms","quasielastic"},string("pair_phi_diff_")+to_string(3)+string("-Bin-") + to_string(bin_num),suffix)
                    .Scale(4).XRange(90,270);
                const hist<> data_copl_mc=
                    Hist(MC,ppn_reaction,{"Histograms","quasielastic"},string("pair_phi_diff_")+to_string(3)+string("-Bin-") + to_string(bin_num),suffix)
                    .Scale(4).XRange(90,270);
                cout << endl << Qmsg<<suffix << endl;
                cout << endl << "Fitting" << endl;

                const auto acc=data_copl_mc.TotalSum()/N;
                const auto res=SystematicError<ppn_fit_range>(
                    [&data_copl,&data_copl_mc,&counter,&Q,&Qmsg,&runmsg,&N,&planarity,&acc]
                    (const double&x){
                    counter++;
                    const auto data_copl_bg=data_copl.XExclude(180.-x,180.+x);
                    Fit<DifferentialMutations<>,ChiSquare,FunctionUncertaintiesEstimation> fit(
                        removeXerorbars(data_copl_bg),
                        [](const ParamSet&X,const ParamSet&P){return Polynom<2>(X[0],P);}
                    );
                    fit.SetUncertaintyCalcDeltas({0.001,0.001,0.001});
                    fit.Init(500,make_shared<InitialDistributions>()
                        <<make_shared<DistribGauss>(0,5000)
                        <<make_shared<DistribGauss>(0,5)
                        <<make_shared<DistribGauss>(0,1)
                    );
                    while(!fit.AbsoluteOptimalityExitCondition(0.000001))fit.Iterate();
                    cout << "Fitting: " << fit.iteration_count() << " iterations; "
                        << fit.Optimality() << "<chi^2<"
                        << fit.Optimality(fit.PopulationSize() - 1)
                        << endl;
                    const auto &P = fit.ParametersWithUncertainties();
                    const auto BG=[&P](const value<>&x){return Polynom<2>(x,P);};
                    const auto subtr=data_copl-BG;
                    const auto summ=subtr.XRange(120,240);
                    const SortedPoints<> simulation_curve=(toLine(data_copl.CloneEmptyBins()+data_copl_mc*summ.TotalSum()/N/acc));
                    if(counter==2){
                        const auto max=to_string(data_copl.TransponateAndSort().right().X().max()*1.5);
                        Plot(Q.Contains(21) ? "ppn-v2-above-data-copl" : (Q.Contains(-39) ? "ppn-v2-below-data-copl" : ""),7)
                            .Hist(data_copl,"DATA").Hist(data_copl_bg)
                            .Hist(data_copl.CloneEmptyBins()+BG,"BG")
                            << "set title '" + runmsg+ Qmsg + "'" <<"set key on"<< "set ylabel 'Events count, n.d.'"
                            << "set yrange [-100:"+max+"]" << "set xlabel " + planarity<< "set xrange [90:270]"<<"set xtics 30";
                        Plot(Q.Contains(21) ? "ppn-v2-above-data-copl-norm" : (Q.Contains(-39) ? "ppn-v2-below-data-copl-norm" : ""),7)
                            .Hist(subtr).Hist(summ,"DATA-BG")
                            .Line(simulation_curve,"Simulation")
                            .Line(Points<>{{subtr.left().X().min(), 0.0},{subtr.right().X().max(), 0.0}})
                            << "set title 'Subtracted background'" <<"set key on"<< "set ylabel 'Events count, n.d.'"
                            << "set yrange [-100:"+max+"]" << "set xlabel " + planarity<< "set xrange [90:270]"<<"set xtics 30";
                    }
                    return extend_value<1,2>(summ.TotalSum())/extend_value<2,2>(acc);
                })();
                return res*trigger_elastic1.scaling/extend_value<2,2>(SIGMA[bin_num].Y());
            };
            luminosity<<make_point(Q,RawSystematicError({pbeam_corr},[&LC,&Q,&lum_bc_m,&lum_bc_p,&lum_bc_z](const string&suffix){
                const auto res=(suffix=="_")?RawSystematicError({ppn_th1,ppn_th2,ppn_t1,ppn_t2},LC)():LC(suffix);
                if(suffix=="_")lum_bc_z<<make_point(Q,res);
                if(suffix=="00-")lum_bc_m<<make_point(Q,res);
                if(suffix=="00+")lum_bc_p<<make_point(Q,res);
                return res;
            })());
        }
    }
    Plot("ppn-v2-acceptance",7)
        .Hist(wrap_hist(acceptance))
            << "set key on" << "set yrange [0:0.2]" << "set xrange [-70:30]"
            << "set xlabel 'Q, MeV/c^2'" << "set ylabel 'Efficiency, n.d.'"<<"set xtics 20";

    Plot("luminosity-v2",3)
        .Hist_2bars<1,2>(lum_bc_z, "","","LUMINOSITYc_z")
        .Hist_2bars<1,2>(lum_bc_p, "","","LUMINOSITYc_p")
        .Hist_2bars<1,2>(lum_bc_m, "","","LUMINOSITYc_m")
            << "set title 'Integrated luminosity " + runmsg + "'"
            << "set key on" << "set xlabel 'Q, MeV/c^2'"<<"set xtics 20"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:]";
    const auto prev_luminosity = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYf"));
    const hist<> sasha=SortedPoints<double,value<>>(Plotter::Instance().GetPoints<double,value<>>("luminosity_khr"));
    const auto hirange=to_string(wrap_hist(luminosity).TransponateAndSort().right().X().max()*1.5);
    Plot("luminosity-v2-compare",3)
        .Hist_2bars<1,2>(luminosity, "ppn_{sp}","","LUMINOSITYc")
        .Hist_2bars<1,2>(prev_luminosity,"3He+eta")
            << "set title '" + runmsg + "'"
            << "set key on" << "set xlabel 'Q, MeV/c^2'"<<"set xtics 20"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:"+hirange+"]";
    Plot("luminosity-v2-compare-light",3)
        .Hist(wrap_hist(luminosity), "ppn_{sp}")
        .Hist(wrap_hist(prev_luminosity),"3He+eta")
            << "set title '" + runmsg + "'"
            << "set key on" << "set xlabel 'Q, MeV/c^2'"<<"set xtics 20"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:"+hirange+"]";
    cout<<"luminosity: "<<luminosity.TotalSum()<<endl;
    Plot("luminosity-v2-compare-estimation",3)
        .Hist(wrap_hist(luminosity)*runs.second/runs.first, "ppn_{sp}")
        .Hist(wrap_hist(prev_luminosity)*runs.second/runs.first, "3He+eta")
        .Hist(sasha, "A. Khreptak")
            << "set key on" << "set xlabel 'Q, MeV/c^2'"<<"set xtics 20"
            << "set ylabel 'Luminosity estimation, nb^{-1}'"
            << "set xrange [-70:30]" << "set yrange [0:100]";
    cout<<"Full luminosity estimation: "<<luminosity.TotalSum()*runs.second/runs.first<<endl;
}

