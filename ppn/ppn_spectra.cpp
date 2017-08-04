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
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_0")/norm.TotalSum().val(),"ppn_{sp}")
    <<"set key on"<<"set title 'Planarity. MC'"<<"set yrange [0:]"<<"set xlabel "+planarity;
    Plot<>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_1")/norm_pd.TotalSum().val(),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_1")/norm.TotalSum().val(),"ppn_{sp}")
    <<"set key on"<<"set title 'Planarity. MC. Cut'"<<"set yrange [0:]"<<"set xlabel "+planarity;
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
    hist<> acceptance,acceptance_pd,chi_sq,chi_sq_mc,luminosity;
    vector<hist<>> fit_params_mc;
    vector<hist<>> fit_params;
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

	const hist<> mc_ppn=Hist(MC,"ppn_qf",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num))
	.Scale(6).XRange(50,250);
	const hist<> mc_pd=Hist(MC,"pd",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num))
	.Scale(6).XRange(50,250);
	const hist<> data=Hist(DATA,"E",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num))
	.Scale(6).XRange(50,250);
	const hist<> nmc_ppn=mc_ppn/N;
	const hist<> nmc_pd=mc_pd/N_pd;
	acceptance<<point<value<>>(Q,nmc_ppn.TotalSum());
	acceptance_pd<<point<value<>>(Q,nmc_pd.TotalSum());
	cout<<endl<<Qmsg<<endl;
	typedef Func4<Novosibirsk,Arg<0>,Par<1>,Par<2>,Par<3>> ToIntegrate;
	typedef Mul<Par<0>,ToIntegrate> Foreground;
	FitFunction<DifferentialMutations<Uncertainty>,Foreground> 
	FitMC(make_shared<FitPoints>(nmc_ppn));
	FitMC.SetFilter([](const ParamSet&P){
	    return (P[0]>0)&&(P[2]>0);
	});
	FitMC.Init(100,make_shared<InitialDistributions>()
	    <<make_shared<DistribUniform>(0.8,1.2)
	    <<make_shared<DistribUniform>(100,120)
	    <<make_shared<DistribUniform>(10,20)
	    <<make_shared<DistribGauss>(0,0.5)
	    ,r_eng
	);
	FitMC.SetUncertaintyCalcDeltas(parEq(Foreground::ParamCount,0.01));
	while(!FitMC.AbsoluteOptimalityExitCondition(0.000001)){
	    FitMC.Iterate(r_eng);
	    cout<<"MC: "<<FitMC.iteration_count()<<" iterations; "
	    <<FitMC.Optimality()<<"<chi^2<"
	    <<FitMC.Optimality(FitMC.PopulationSize()-1)
	    <<"          \r";
	}
	const auto&P0=FitMC.ParametersWithUncertainties();
	for(size_t i=0;i<Foreground::ParamCount;i++){
	    if(fit_params_mc.size()==i)fit_params_mc.push_back(hist<>());
	    fit_params_mc[i]<<point<value<>>(Q,P0[i]);
	}
	cout<<endl;
	Plot<>().Hist(nmc_ppn,"MC")
	.Line(SortedPoints<>([&FitMC](double x){return FitMC({x});},ChainWithStep(50.,1.,250.)),"fit")
	<<"set key on"<<"set title 'MC "+Qmsg+"'"<<"set yrange [0:]"
	<<"set xlabel "+thth<<"set ylabel 'counts normalized'";
	static const size_t fgp=Foreground::ParamCount;
	typedef Mul<
	    Func3<FermiFunc,Arg<0>,Par<fgp+0>,Par<fgp+1>>,
	    PolynomFunc<Arg<0>,fgp+2,1>
	>BackGround;
	const auto&data_count=data.TotalSum().val();
	FitFunction<DifferentialMutations<Uncertainty>,Add<Foreground,BackGround>> 
	FitData(make_shared<FitPoints>(data));
	FitData.SetFilter([](const ParamSet&P){
	    return (P[0]>0)&&(P[1]>100)&&(P[1]<120)&&(P[2]>0)
	    &&(P[fgp+0]>50)&&(P[fgp+0]<80)&&(P[fgp+1]<0);
	});
	FitData.Init(200,make_shared<InitialDistributions>()
	    <<make_shared<DistribUniform>(0,data_count)
	    <<make_shared<DistribUniform>(100,120)
	    <<make_shared<DistribUniform>(10,20)
	    <<make_shared<DistribGauss>(0,0.5)
	    <<make_shared<DistribUniform>(60,70)
	    <<make_shared<DistribUniform>(-5,0)
	    <<make_shared<DistribUniform>(0,0.01*data_count)
	    <<make_shared<DistribUniform>(-100,0)
	    ,r_eng
	);
	FitData.SetUncertaintyCalcDeltas(parEq(BackGround::ParamCount,0.01));
	while(!FitData.AbsoluteOptimalityExitCondition(0.000001)){
	    FitData.Iterate(r_eng);
	    cout<<"DATA: "<<FitData.iteration_count()<<" iterations; "
	    <<FitData.Optimality()<<"<chi^2<"
	    <<FitData.Optimality(FitData.PopulationSize()-1)
	    <<"          \r";
	}
	const auto&P1=FitData.ParametersWithUncertainties();
	const auto&p1=FitData.Parameters();
	for(size_t i=0;i<BackGround::ParamCount;i++){
	    if(fit_params.size()==i)fit_params.push_back(hist<>());
	    fit_params[i]<<point<value<>>(Q,P1[i]);
	}
	Plot<>().Hist(data,"DATA")
	.Line(SortedPoints<>([&FitData](double x){return FitData({x});},ChainWithStep(50.,1.,250.)),"fit")
	.Line(SortedPoints<>([&p1](double x){return BackGround()({x},p1);},ChainWithStep(50.,1.,250.)),"bg")
	<<"set key on"<<"set title 'Data "+Qmsg+"'"<<"set yrange [0:]"
	<<"set xlabel "+thth<<"set ylabel 'counts'";
	
	chi_sq_mc<<point<value<>>(Q,FitMC.Optimality()/(nmc_ppn.size()-FitMC.ParamCount()));
	chi_sq<<point<value<>>(Q,FitData.Optimality()/(data.size()-FitData.ParamCount()));
	const double 
	I0=Sympson([&FitMC](const double&x){return ToIntegrate()({x},FitMC.Parameters());},50.,250.,0.001),
	I1=Sympson([&FitData](const double&x){return ToIntegrate()({x},FitData.Parameters());},50.,250.,0.001);
	luminosity << point<value<>>(Q,
	    ((P1[0]*I1)/(P0[0]*I0)/SIGMA[bin_num].Y())*double(trigger_elastic1.scaling)
	);
    }
    Plot<>().Hist(acceptance,"ppn_{sp}").Hist(acceptance_pd,"pd")<<"set key on"
    <<"set title 'Acceptance'"<<"set yrange [0:]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
    for(size_t i=0;i<fit_params_mc.size();i++){
	Plot<double>().Hist(fit_params_mc[i])
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'mc parameter"+to_string(i)+"'";
    }
    for(size_t i=0;i<fit_params.size();i++){
	Plot<double>().Hist(fit_params[i])
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'parameter"+to_string(i)+"'";
    }
    Plot<>().Hist(fit_params[1],"DATA").Hist(fit_params_mc[1],"MC")
    << "set xlabel 'Q, MeV'"<< "set ylabel 'Position'"<<"set key on"
    <<"set yrange [50:150]";
    Plot<>().Hist(fit_params[2],"DATA").Hist(fit_params_mc[2],"MC")
    << "set xlabel 'Q, MeV'"<< "set ylabel 'Width'"<<"set key on"
    <<"set yrange [0:30]";
    Plot<>().Hist(fit_params[3],"DATA").Hist(fit_params_mc[3],"MC")
    << "set xlabel 'Q, MeV'"<< "set ylabel 'Asymmetry'"<<"set key on"
    <<"set yrange [-0.2:1.0]";
    
    Plot<double>().Hist(chi_sq_mc,"MC").Hist(chi_sq,"DATA")
    << "set xlabel 'Q, MeV'" <<"set key on"
    << "set ylabel 'chi^2/d, n.d.'" 
    << "set yrange [0:]"<<"unset log y";

    Plot<double>().Hist(luminosity) 
    << "set title 'Integrated luminosity estimation ("+runmsg+")'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integrated luminosity, nb^{-1}'" 
    << "set xrange [-70:30]"<< "set yrange [0:]";
}
