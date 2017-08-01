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
#include <Genetic/searchmin.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
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
const BiSortedPoints<> ReadCrossSection(){
    BiSortedPoints<> result(ChainWithCount(181,0.,PI<>()),ChainWithCount(13,1.000,2.200));
    for(size_t degree=0;degree<=180;degree++){
	ifstream file("crosssections/Theta_"+to_string(degree)+".txt");
	for(double E=0,C=0;(file>>E>>C);result.Bin(degree,(size_t(E)-1000)/100)=C);
	file.close();
    }
    return result;
}
const SortedPoints<> IntegrateCrossSection(const BiSortedPoints<>&angular){
    SortedPoints<> result;
    for(size_t index=0;index<angular.Y().size();index++){
	const auto&E=angular.Y()[index];
	const auto IntT=Int_Trapez_Table(angular.CutX(index)*[](const double&th){return sin(th);});
	result<<point<>(E,IntT[IntT.size()-1].Y());
    }
    return result;
}
const SortedPoints<value<>> ConvertCrossSections(const SortedPoints<>&momentum){
    LinearInterpolation<value<>> q_form;
    static const Reaction he3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
    for(const auto&P:momentum){
	q_form<<point<value<>>(he3eta.P2Q(P.X())*1000.,P.Y()*1000000);
    }
    return SortedPoints<value<>>(q_form.func(),BinsByStep(-70.,2.5,+30.));
}
int main(){
    const auto runs=PresentRuns("E");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";
    const string th1="'Theta_1, deg'",th2="'Theta_2, deg'",e1="'Edep_1, GeV'",e2="'Edep_2, GeV'",
    thth="'Theta_1+1.6Theta_2, deg'",planarity="'|Phi_1-Phi_2-180^o|, deg'";
    const hist<> norm=Hist(MC,"ppn_qf",{"Histograms","elastic"},"0-Reference");
    const hist<> norm_pd=Hist(MC,"pd",{"Histograms","elastic"},"0-Reference");
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"ppn");
    Plot<>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_0")/norm_pd.TotalSum().val(),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_0")/norm.TotalSum().val(),"ppn_{sp}")<<"set key on"
    <<"set title 'Planarity. MC'"<<"set yrange [0:]"<<"set xlabel "+planarity;
    Plot<>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_1")/norm_pd.TotalSum().val(),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_1")/norm.TotalSum().val(),"ppn_{sp}")<<"set key on"
    <<"set title 'Planarity. MC. Cut'"<<"set yrange [0:]"<<"set xlabel "+planarity;
    Plot<>()
    .Hist(Hist(DATA,"E",{"Histograms","elastic"},"pair_phi_diff_0"))
    .Hist(Hist(DATA,"E",{"Histograms","elastic"},"pair_phi_diff_1"))
    <<"set title 'Planarity. Data "+runmsg+"'"<<"set yrange [0:]"<<"set xlabel "+planarity;

    PlotHist2d<>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'MC pd'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'MC ppn_{sp}'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA,"E",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'Data "+runmsg+"'"<<"set xlabel "+th1<<"set ylabel "+th2;


    PlotHist2d<>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC pd'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC ppn_{sp}'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA,"E",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. Data "+runmsg+"'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC pd'"<<"set xlabel "+th1<<"set ylabel "+e1;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1 MC ppn_{sp}'"<<"set xlabel "+th1<<"set ylabel "+e1;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA,"E",{"Histograms","elastic"},"t_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. Data "+runmsg+"'"<<"set xlabel "+th1<<"set ylabel "+e1;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"e_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC pd'"<<"set xlabel "+e1<<"set ylabel "+e2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"e_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC ppn_{sp}'"<<"set xlabel "+e1<<"set ylabel "+e2;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA,"E",{"Histograms","elastic"},"e_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. Data "+runmsg+"'"<<"set xlabel "+e1<<"set ylabel "+e2;
    Plot<>()
    .Line(Hist(MC,"pd",{"Histograms","elastic"},"theta_sum_22-AllBins").toLine()/norm_pd.TotalSum().val(),"pd")
    .Line(Hist(MC,"ppn_qf",{"Histograms","elastic"},"theta_sum_22-AllBins").toLine()/norm.TotalSum().val(),"ppn_{sp}")
    <<"set title 'MC'"<<"set key on"<<"set yrange [0:]"<<"set xlabel "+thth;
    Plot<>()
    .Hist(Hist(DATA,"E",{"Histograms","elastic"},"theta_sum_22-AllBins"))
    <<"set title 'Data "+runmsg+"'"<<"set key on"<<"set yrange [0:]"<<"set xlabel "+thth;


    RANDOM r_eng;
    hist<> acceptance,acceptance_pd,chi_sq,luminosity,luminosity2;
    vector<hist<>> fit_params{hist<>(),hist<>(),hist<>(),hist<>()};
    const auto diff_cs=ReadCrossSection();
    PlotHist2d<>(sp2).Surface(diff_cs.Clone().FullCycleVar([](double&z){z=log10(z);}));
    const auto p_cs=IntegrateCrossSection(diff_cs);
    Plot<>().Line(p_cs);
    const auto SIGMA=ConvertCrossSections(p_cs);
    Plot<double>().Hist(SIGMA)
    << "set title 'ppn_{sp} cross section'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'cross section, nb'" 
    << "set xrange [-70:30]"<< "set yrange [0:]";
    for(size_t bin_num=0,bin_count=norm.size();bin_num<bin_count;bin_num++){
	const auto&Q=norm[bin_num].X();
	const auto&N=norm[bin_num].Y();
	const auto&N_pd=norm_pd[bin_num].Y();
	const string Qmsg=static_cast<stringstream&>(stringstream()
	    <<"Q in ["<<setprecision(3)
	    <<Q.min()<<"; "<<Q.max()<<"] MeV"
	).str();

	const hist<> mc_ppn=Hist(MC,"ppn_qf",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num)).Scale(4).XRange(80,225);
	const hist<> mc_pd=Hist(MC,"pd",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num)).Scale(4).XRange(80,225);
	const hist<> nmc_ppn=mc_ppn/N;
	const hist<> nmc_pd=mc_pd/N_pd;
	acceptance<<point<value<>>(Q,mc_ppn.TotalSum()/N);
	acceptance_pd<<point<value<>>(Q,mc_pd.TotalSum()/N_pd);
	const hist<> data=Hist(DATA,"E",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num)).Scale(4).XRange(80,225);
	Plot<>().Hist(nmc_ppn/(mc_ppn[0].X().uncertainty()*2.),"ppn_{sp}")
	.Hist(nmc_pd/(mc_pd[0].X().uncertainty()*2.),"pd")
	<<"set key on"<<"set title 'MC "+Qmsg+"'"<<"set yrange [0:0.015]"
	<<"set xlabel "+thth<<"set ylabel 'Differential acceptance, deg^{-1}'";
	
	SearchMin<DifferentialMutations<Uncertainty>>
	fit([&nmc_pd,&nmc_ppn,&data](const ParamSet&P){
	    double res=0;
	    for(size_t i=1,n=data.size();i<n;i++){
		const value<double> exp_p=data[i].Y(),
		the_p=nmc_ppn[i].Y()*P[0]+nmc_pd[i].Y()*P[1]+P[2]+data[i].X().val()*P[3];
		res+=exp_p.NumCompare(the_p);
	    }
	    return res;
	});
	fit.SetMutationCoefficient(0.8);
	fit.SetFilter([](const ParamSet&P)->bool{return (P[0]>0)&&(P[1]>0);});
	const auto&data_count=data.TotalSum().val();
	fit.Init(100,
	    make_shared<InitialDistributions>()
		<<make_shared<DistribUniform>(0.0,2.0*data_count)
		<<make_shared<DistribUniform>(0.0,2.0*data_count)
		<<make_shared<DistribGauss>(0.0,0.02*data_count)
		<<make_shared<DistribGauss>(0.0,0.02*data_count)
	    ,r_eng
	);
	while(!fit.AbsoluteOptimalityExitCondition(0.000000000001)){
	    fit.Iterate(r_eng);
	    cout<<fit.iteration_count()<<" iterations; "
	    <<fit.Optimality()<<"<chi^2<"
	    <<fit.Optimality(fit.PopulationSize()-1)
	    <<"          \r";
	}
	fit.SetUncertaintyCalcDeltas({0.1,0.1,0.1,0.1});
	const auto&P=fit.ParametersWithUncertainties();
	for(size_t i=0;i<fit_params.size();i++){
	    fit_params[i]<< point<value<>>(Q,P[i]);
	}
	chi_sq << point<value<>>(Q,fit.Optimality()/(data.size()-fit.ParamCount()));
	const hist<> 
	PPN=nmc_ppn*P[0],PD=nmc_pd*P[1],
	BG=data.Clone().Transform([&P](const value<>&x,const value<>&){return P[2]+P[3]*x.val();});
	Plot<>().Hist(data,"data").Line(hist<>(PPN+PD+BG).toLine(),"ppn+pd+bg")
	.Line(BG.toLine(),"bg").Line(hist<>(PPN+BG).toLine(),"ppn+bg")
	<<"set key on"<<"set title 'Data "+Qmsg+" "+runmsg+"'"<<"set yrange [0:]"
	<<"set xlabel "+thth<<"set ylabel 'Count'";
	luminosity << point<value<>>(Q,(P[0]/SIGMA[bin_num].Y())*double(trigger_elastic.scaling));
	const hist<> FG=data-BG-PD;
	Plot<>().Hist(FG)
	<<"set key on"<<"set title 'Foreground estimation "+Qmsg+" "+runmsg+"'"<<"set yrange [0:]"
	<<"set xlabel "+thth<<"set ylabel 'Count'";
	luminosity2 << point<value<>>(Q,(FG.TotalSum()/nmc_ppn.TotalSum()/SIGMA[bin_num].Y()) * trigger_elastic.scaling);
    }
    Plot<>().Hist(acceptance,"ppn_{sp}").Hist(acceptance_pd,"pd")<<"set key on"
    <<"set title 'Acceptance'"<<"set yrange [0:]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
    for(size_t i=0;i<fit_params.size();i++)
	Plot<double>().Hist(fit_params[i])
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'parameter"+to_string(i)+"'";
    Plot<double>().Hist(chi_sq)
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'chi^2/d, n.d.'" 
    << "set yrange [0:]"<<"unset log y";

    Plot<double>().Hist(luminosity) 
    << "set title 'Integrated luminosity estimation ("+runmsg+")'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integrated luminosity, nb^{-1}'" 
    << "set xrange [-70:30]"<< "set yrange [0:]";
    Plot<double>().Hist(luminosity2) 
    << "set title 'Integrated luminosity estimation ("+runmsg+")'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integrated luminosity, nb^{-1}'" 
    << "set xrange [-70:30]"<< "set yrange [0:]";
}
