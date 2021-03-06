// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <math_h/sigma3.h>
#include <math_h/lorentzvector.h>
#include <math_h/interpolate.h>
#include <Genetic/searchmin.h>
#include <Genetic/uncertainties.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Parameters/parameters.h>
#include <Parameters/systematic.h>
#include <Kinematics/particles.h>
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
const value<> branching_ratio1{0.3931,0.0020};
const value<> branching_ratio2{0.3257,0.0023};
const list<size_t> params{
    pbeam_corr,he3_cut_h,he3_theta_cut,
    gamma_E_thr,time_dt,time_t1,time_t2,eta_theta_thr,he3mm_cut,
    gamma_mm_lo,gamma_mm_hi,gamma_im_lo,gamma_im_hi,
    gamma_im_lo6,gamma_im_hi6,three_pi0
};
const vector<string> histpath_central_reconstr1 = {"Histograms", "He3nCentralGammas2"};
const vector<string> reaction1 = {"bound1-2g","bound2-2g","bound3-2g"};
const vector<string> histpath_central_reconstr2 = {"Histograms", "He3nCentralGammas6"};
const vector<string> reaction2 = {"bound1-6g","bound2-6g","bound3-6g"};
class LuminosityGetter{
private:
    ext_hist<2> khreptak_data;
public:
    LuminosityGetter(){
        cout<<"getting luminosity"<<endl;
	const auto data=Plotter::Instance().GetPoints<>("khreptak-luminosity");
	value<> q(-70.0+1.25,1.25);
	for(const auto&p:data){
		khreptak_data<<make_point(q,uncertainties(p.Y(),0.0,0.0));
		q+=2.5;
	}
    }
    ~LuminosityGetter(){}
    static const LuminosityGetter&Instance(){
        static LuminosityGetter instance;
        return instance;
    }
    const ext_hist<2>&get(const string&suffix)const{
        return khreptak_data;
    }
};
typedef SearchMin<DifferentialMutations<>,UncertaintiesEstimation> Fitter;
Cache<string,ext_hist<2>> DATA1,DATA2;
const auto Qbins=BinsByStep(-67.5,2.5,0.0);
Fitter BGfit(const string&suffix,bool plot=false){
    const auto&lum=LuminosityGetter::Instance().get(suffix);
    cout<<"Histogram for "<<suffix<<": "<<"data1"<<endl;
    const auto&data1=DATA1(suffix,[&suffix,&lum](){
        return ext_hist<2>([&suffix,&lum](const value<>&Q){
                const size_t bin_num=trunc((Q.val()+70.)/2.5);
                cout<<"Preparing histogram for "<<suffix<<": "<<"data1["<<Q.min()<<":"<<Q.max()<<"]MeV"<<endl;
                const auto ac=SystematicError<bound_state_reaction_index>([&bin_num,&suffix](const int i){
                    const auto &r = reaction1[i];
                    const auto MC_TIM=Hist(MC, r, histpath_central_reconstr1, "TIM6-Bin-"+to_string(bin_num),suffix).XRange(-0.3,0.3);
                    static Cache<string,hist<>> NORM;
                    const auto &N = NORM(r+suffix,[&r,&suffix](){return Hist(MC, r, histpath_central_reconstr1, "0-Reference",suffix);})[bin_num].Y();
                    return extend_value<2,2>(std_error(MC_TIM.TotalSum().val())/N);
                })();
            const auto TIM=Hist(DATA, "All", histpath_central_reconstr1, string("TIM6-Bin-") + to_string(bin_num),suffix).XRange(-0.3,0.3);
            return (extend_value<1,2>(std_error(TIM.TotalSum().val()))*trigger_he3_forward.scaling)/(ac*lum[bin_num].Y());
        },Qbins);
    });
    cout<<"Histogram for "<<suffix<<": "<<"data2"<<endl;
    const auto&data2=DATA2(suffix,[&suffix,&lum](){
        return ext_hist<2>([&suffix,&lum](const value<>&Q){
            const size_t bin_num=trunc((Q.val()+70.)/2.5);
            cout<<"Preparing histogram for "<<suffix<<": "<<"data2["<<Q.min()<<":"<<Q.max()<<"]MeV"<<endl;
            const auto ac=SystematicError<bound_state_reaction_index>([&bin_num,&suffix](const int i){
                const auto &r = reaction2[i];
                const auto MC_TIM=Hist(MC, r, histpath_central_reconstr2, "TIM7-Bin-"+to_string(bin_num),suffix).XRange(-0.3,0.3);
                static Cache<string,hist<>> NORM;
                const auto &N = NORM(r+suffix,[&r,&suffix](){return Hist(MC, r, histpath_central_reconstr2, "0-Reference",suffix);})[bin_num].Y();
                return extend_value<2,2>(std_error(MC_TIM.TotalSum().val())/N);
            })();
            const auto TIM=Hist(DATA, "All", histpath_central_reconstr2, string("TIM7-Bin-") + to_string(bin_num),suffix).XRange(-0.3,0.3);
            return (extend_value<1,2>(std_error(TIM.TotalSum().val()))*trigger_he3_forward.scaling)/(ac*lum[bin_num].Y());
        },Qbins);
    });
    const auto Data1=take_uncertainty_component<1>(data1),
               Data2=take_uncertainty_component<1>(data2);
    cout<<"Fitting for BG"<<endl;
    Fitter FIT([Data1,Data2](const ParamSet&P){
        const auto F1=[&P](const double&x){return x*P[0]+P[1];};
        const auto F2=[&P](const double&x){return x*P[2]+P[3];};
        double res=0;
        for(size_t i=0;i<Qbins.size();i++){
            res+=Data1[i].Y().NumCompare(F1(Qbins[i].val()));
            res+=Data2[i].Y().NumCompare(F2(Qbins[i].val()));
        }
        return res;
    });
    FIT.SetFilter([](const ParamSet&P){
	return (pow(P[0],2)<1)&&(pow(P[2],2)<1)&&
	    (P[1]>0)&&(P[1]<50)&&(P[3]>0)&&(P[3]<50);
    });
    FIT.Init(100,make_shared<InitialDistributions>()
        <<make_shared<DistribUniform>(-0.1,0)
        <<make_shared<DistribUniform>(10,30)
        <<make_shared<DistribUniform>(-0.1,0)
        <<make_shared<DistribUniform>(20,40)
    );
    while(!FIT.AbsoluteOptimalityExitCondition(0.0000001))FIT.Iterate();
    cout << "Fit result: " << FIT.iteration_count() << " iterations; "
        << FIT.Optimality() << "<chi^2<"
        << FIT.Optimality(FIT.PopulationSize() - 1)
        << endl;
    FIT.SetUncertaintyCalcDeltas({0.001,0.1,0.001,0.1});
    const auto&P=FIT.Parameters();
    if(plot){
        const auto background1=toLine(hist<>([&P](const value<>&x){return x*P[0]+P[1];},Qbins));
        const auto background2=toLine(hist<>([&P](const value<>&x){return x*P[2]+P[3];},Qbins));
        Plot("UpperLimit2-LinearFit1",5)
            .Hist(Data1).Line(background1)
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:60]"
                << "set title 'pd->^3He2γ'"<<"set key right top";
        Plot("UpperLimit2-LinearFit2",5)
            .Hist(Data2).Line(background2)
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:60]"
                << "set title 'pd->^3He6γ'"<<"set key right top";
    }
    return FIT;

}
Fitter BWfit(const string&suffix,size_t power,const double&B,const double&G,bool plot=false){
    const auto&lum=LuminosityGetter::Instance().get(suffix);
    cout<<"Histogram for "<<suffix<<": "<<"data1"<<endl;
    const auto&data1=DATA1(suffix,[&suffix,&lum](){
        return ext_hist<2>([&suffix,&lum](const value<>&Q){
                const size_t bin_num=trunc((Q.val()+70.)/2.5);
                cout<<"Preparing histogram for "<<suffix<<": "<<"data1["<<Q.min()<<":"<<Q.max()<<"]MeV"<<endl;
                const auto ac=SystematicError<bound_state_reaction_index>([&bin_num,&suffix](const int i){
                    const auto &r = reaction1[i];
                    const auto MC_TIM=Hist(MC, r, histpath_central_reconstr1, "TIM6-Bin-"+to_string(bin_num),suffix).XRange(-0.3,0.3);
                    static Cache<string,hist<>> NORM;
                    const auto &N = NORM(r+suffix,[&r,&suffix](){return Hist(MC, r, histpath_central_reconstr1, "0-Reference",suffix);})[bin_num].Y();
                    return extend_value<2,2>(std_error(MC_TIM.TotalSum().val())/N);
                })();
            const auto TIM=Hist(DATA, "All", histpath_central_reconstr1, string("TIM6-Bin-") + to_string(bin_num),suffix).XRange(-0.3,0.3);
            return (extend_value<1,2>(std_error(TIM.TotalSum().val()))*trigger_he3_forward.scaling)/(ac*lum[bin_num].Y());
        },Qbins);
    });
    cout<<"Histogram for "<<suffix<<": "<<"data2"<<endl;
    const auto&data2=DATA2(suffix,[&suffix,&lum](){
        return ext_hist<2>([&suffix,&lum](const value<>&Q){
            const size_t bin_num=trunc((Q.val()+70.)/2.5);
            cout<<"Preparing histogram for "<<suffix<<": "<<"data2["<<Q.min()<<":"<<Q.max()<<"]MeV"<<endl;
            const auto ac=SystematicError<bound_state_reaction_index>([&bin_num,&suffix](const int i){
                const auto &r = reaction2[i];
                const auto MC_TIM=Hist(MC, r, histpath_central_reconstr2, "TIM7-Bin-"+to_string(bin_num),suffix).XRange(-0.3,0.3);
                static Cache<string,hist<>> NORM;
                const auto &N = NORM(r+suffix,[&r,&suffix](){return Hist(MC, r, histpath_central_reconstr2, "0-Reference",suffix);})[bin_num].Y();
                return extend_value<2,2>(std_error(MC_TIM.TotalSum().val())/N);
            })();
            const auto TIM=Hist(DATA, "All", histpath_central_reconstr2, string("TIM7-Bin-") + to_string(bin_num),suffix).XRange(-0.3,0.3);
            return (extend_value<1,2>(std_error(TIM.TotalSum().val()))*trigger_he3_forward.scaling)/(ac*lum[bin_num].Y());
        },Qbins);
    });
    const auto Data1=take_uncertainty_component<1>(data1),
               Data2=take_uncertainty_component<1>(data2);
    cout<<"Fitting for "<<suffix<<": B="<<B<<"; G="<<G<<endl;
    const auto static Q2P=[](const double&Q){
	const auto TM=Q*0.001+Particle::eta().mass()+Particle::he3().mass();
	const auto cm=binaryDecay(TM,Particle::p().mass(),Particle::d().mass());
	const auto p=cm.first.Transform(cm.second.Beta());
	return p.P().M();
    };
    const hist<> BW([&B,&G](const value<>&Q)->value<>{return BreitWigner(Q.val(),B,G)/BreitWigner(B,B,G);},Qbins);
    Fitter FIT([BW,Data1,Data2,power](const ParamSet&P){
        auto F1=BW*branching_ratio1*P[0]+[&P](const value<>&x)->value<>{return x*P[1]+P[2];};
        auto F2=BW*branching_ratio2*P[0]+[&P](const value<>&x)->value<>{return x*P[3]+P[4];};
	if(power==2){
	    F1+=[&P](const value<>&x)->value<>{return x*x*P[5];};
	    F2+=[&P](const value<>&x)->value<>{return x*x*P[6];};
	}
        double res=0;
        for(size_t i=0;i<Qbins.size();i++){
            res+=Data1[i].Y().NumCompare(F1[i].Y());
            res+=Data2[i].Y().NumCompare(F2[i].Y());
        }
        return res;
    });
    FIT.SetFilter([](const ParamSet&P){
	    return (P[0]>=0)&&
		(pow(P[1],2)<1)&&(pow(P[3],2)<1)&&
		(P[2]>0)&&(P[2]<50)&&(P[4]>0)&&(P[4]<50);
    });
    auto init=make_shared<InitialDistributions>()
        <<make_shared<DistribUniform>(0,10)
        <<make_shared<DistribUniform>(-0.1,0)
        <<make_shared<DistribUniform>(10,30)
        <<make_shared<DistribUniform>(-0.1,0)
        <<make_shared<DistribUniform>(20,40);
    if(power==2)init<<make_shared<DistribUniform>(-0.01,0.01)<<make_shared<DistribUniform>(-0.01,0.01);
    FIT.Init(100,init);
    while(!FIT.AbsoluteOptimalityExitCondition(0.0000001))FIT.Iterate();
    cout << "Fit result: " << FIT.iteration_count() << " iterations; "
        << FIT.Optimality() << "<chi^2<"
        << FIT.Optimality(FIT.PopulationSize() - 1)
        << endl;
    FIT.SetUncertaintyCalcDeltas({0.1,0.001,0.1,0.001,0.1});
    const auto&P=FIT.Parameters();
    if(plot){
        const auto background1=toLine(hist<>([&P](const value<>&x){return x*P[1]+P[2];},Qbins));
        const auto background2=toLine(hist<>([&P](const value<>&x){return x*P[3]+P[4];},Qbins));
        const auto fit1=toLine(BW)*branching_ratio1.val()*P[0]+background1;
        const auto fit2=toLine(BW)*branching_ratio2.val()*P[0]+background2;
        const string msg="B="+to_string(B)+"; Г="+to_string(G);
        Plot("UpperLimit2-"+to_string(int(B*10))+"-"+to_string(int(G*10))+"Fit1",5)
            .Hist(Data1).Line(fit1).Line(background1)
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:60]"
                << "set title 'pd->^3He2γ "+msg+"'"<<"set key right top";
        Plot("UpperLimit2-"+to_string(int(B*10))+"-"+to_string(int(G*10))+"Fit2",5)
            .Hist(Data2).Line(fit2).Line(background2)
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:60]"
                << "set title 'pd->^3He6γ "+msg+"'"<<"set key right top";
    }
    return FIT;
}
const double K=1.64485;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "estimation2");
    const auto binding=BinsByStep(-65.0,2.5,0.0),gamma=BinsByStep(0.0,2.5,40.0);
    BiSortedPoints<value<>,value<>,Fitter> fit_results(binding,gamma,[](const ParamSet&){return 0.0;});
    cout<<"Fitting"<<endl;
    fit_results.FullCycleVar([](const value<>&B,const value<>&G,Fitter&res){
	res=BWfit("_",1,B.val(),G.val(),G.Contains(29)||G.Contains(19)||G.Contains(9));
    });
    BiSortedPoints<value<>,value<>,value<>> chisquare(binding,gamma),upper(binding,gamma);
    cout<<"converting"<<endl;
    Plot chisq_plot("UpperLimit2-chisq-1d");
    chisq_plot<<"set xlabel 'Peak position, MeV'"<<"set yrange [0.5:1.3]"
	    <<"set ylabel 'χ^2/d, n.d.'"<<"set key on";
    for(size_t g=0;g<gamma.size();g++){
	SortedPoints<> chi_sq_curve;
	for(size_t b=0;b<binding.size();b++){
	    const auto&src=fit_results[b][g];
	    const auto chisq=src.Optimality()/(Qbins.size()+Qbins.size()-src.ParamCount());
    	    chisquare.Bin(b,g)=chisq;
	    chi_sq_curve<<make_point(binding[b].val(),chisq);
	    const auto&A=src.ParametersWithUncertainties()[0];
	    upper    .Bin(b,g)=A.uncertainty()*K;
	}
	if(g%4==3)chisq_plot.Line(chi_sq_curve,"Width="+to_string(gamma[g].val())+"MeV");
    }
    cout<<"plotting"<<endl;
    PlotHist2d(sp2,"UpperLimit2-chisq").Distr(chisquare)<<"set colorbox"
	<<"set xlabel 'Peak position, MeV'"<<"set ylabel 'Width, MeV'"
	<<"set title 'χ^2/d'";
    PlotHist2d(sp2,"UpperLimit2-upperlimit").Distr(upper)<<"set colorbox"
	<<"set xlabel 'Peak position, MeV'"<<"set ylabel 'Width, MeV'"
	<<"set title 'Upper limit, nb'";
    for(const auto&G:gamma){
        ext_hist<2> A_hist;
        SortedPoints<> upper_limit,systematics;
        for(const auto&B:binding){
	    RawSystematicError calc(params,[&B,&G](const string&suffix){
		return SystematicError<upper_limit_fit_power>([&B,&G,&suffix](const double&power){
		    const auto fit=BWfit(suffix,size_t(power),B.val(),G.val(),false);
		    const auto&A=fit.ParametersWithUncertainties()[0];
        	    return uncertainties(A.val(),A.uncertainty()*K,0.0);
		})();
	    });
	    const auto A=calc();
	    A_hist<<make_point(B.val(),A);
	    upper_limit<<make_point(B.val(),A.uncertainty<1>());
	    systematics<<make_point(B.val(),A.uncertainty<2>());
        }
	Plot("UpperLimit2-Gconst"+to_string(int(G.val()*10)))
	    .Line(upper_limit,"Upper limit")
	    .Line(systematics,"Systematics")
	    <<"set title 'Width = "+to_string(G.val())+" MeV'"<<"set key on"
	    <<"set xlabel 'Peak position, MeV'"<<"set yrange [0:30]"
	    <<"set ylabel 'σ, nb'";
	Plot("UpperLimit2-Gconst"+to_string(int(G.val()*10))+"-Parameter")
	    .Hist(take_uncertainty_component<1>(A_hist),"Fit")
	    .Line(systematics,"Systematics")
	    <<"set title 'Width = "+to_string(G.val())+" MeV'"<<"set key on"
	    <<"set xlabel 'Peak position, MeV'"<<"set yrange [-1:30]"
	    <<"set ylabel 'σ, nb'";
    }
    {
	RawSystematicError calc(params,[](const string&suffix){
	    const auto fit=BGfit(suffix,suffix=="_");
            return uncertainties(fit.Optimality()/(Qbins.size()+Qbins.size()-fit.ParamCount()),0.0,0.0);
        });
	const auto chisq=wrap_value(calc());
	chisq_plot.Line(Points<>{{-70,chisq.val()},{0,chisq.val()}},"Linear fit");
	const hist<> chisq_hist=hist<>()<<make_point(0.0,chisq);
	Plot("UpperLimit2-LinearFitChisq")
	    .Hist(chisq_hist,"","LinearFitChisq")
	    <<"set ylabel 'χ^2/d, n.d.'"
	    <<"set xlabel 'Parameter index'";
    }
}