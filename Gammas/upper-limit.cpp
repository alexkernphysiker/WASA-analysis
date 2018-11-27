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
    ext_hist<2> normal,upper,lower;
public:
    LuminosityGetter(){
        cout<<"getting luminosity"<<endl;
        normal = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc_z"));
        upper = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc_p"));
        lower = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc_m"));
    }
    ~LuminosityGetter(){}
    static const LuminosityGetter&Instance(){
        static LuminosityGetter instance;
        return instance;
    }
    const ext_hist<2>&get(const string&suffix)const{
        return (suffix=="00+")?upper:(suffix=="00-")?lower:normal;
    }
};
typedef SearchMin<DifferentialMutations<>,UncertaintiesEstimation> Fitter;
Cache<string,ext_hist<2>> DATA1,DATA2;
const auto Qbins=BinsByStep(-70.0,2.5,2.5);
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
    const auto Data1=take_uncertainty_component<1>(data1.XRange(-67.5,0.0)),
               Data2=take_uncertainty_component<1>(data2.XRange(-67.5,0.0));
    cout<<"Fitting for BG"<<endl;
    Fitter FIT([Data1,Data2](const ParamSet&P){
        const auto F1=[&P](const double&x){return x*P[0]+P[1];};
        const auto F2=[&P](const double&x){return x*P[2]+P[3];};
        double res=0;
        for(size_t i=0;i<Data1.size();i++)res+=Data1[i].Y().NumCompare(F1(Data1[i].X().val()));
        for(size_t i=0;i<Data2.size();i++)res+=Data2[i].Y().NumCompare(F2(Data2[i].X().val()));
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
        Plot("UpperLimit-LinearFit1",5)
            .Hist(Data1,"data").Line(background1,"linear","","lc rgb \"black\"")
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [10:40]"
                << "set title 'pd->^3He2γ'"<<"set key right top";
        Plot("UpperLimit-LinearFit2",5)
            .Hist(Data2,"data").Line(background2,"linear","","lc rgb \"black\"")
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [15:50]"
                << "set title 'pd->^3He6γ'"<<"set key right top";
    }
    return FIT;

}
Fitter BWfit(const string&suffix,size_t power,
			const double&la,const double&lb,
			const double&B,const double&G,bool plot=false
){
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
    const auto Data1=take_uncertainty_component<1>(data1.XRange(-67.5+la,0.0+lb)),
               Data2=take_uncertainty_component<1>(data2.XRange(Data1.left().X().min(),Data1.right().X().max()));
    cout<<"Fitting for "<<suffix<<": B="<<B<<"; G="<<G<<";A="<<la<<";B="<<lb<<endl;
    const hist<> BW([&B,&G](const value<>&Q)->value<>{return BreitWigner(Q.val(),B,G)/BreitWigner(B,B,G);},Data1);
    Fitter FIT([BW,Data1,Data2,power](const ParamSet&P){
        auto F1=BW*branching_ratio1*P[0]+[&P](const value<>&x)->value<>{return x*P[1]+P[2];};
        auto F2=BW*branching_ratio2*P[0]+[&P](const value<>&x)->value<>{return x*P[3]+P[4];};
	if(power==2){
	    F1+=[&P](const value<>&x)->value<>{return x*x*P[5];};
	    F2+=[&P](const value<>&x)->value<>{return x*x*P[6];};
	}
        double res=0;
        for(size_t i=0;i<Data1.size();i++){
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
        hist<> F1([&P](const value<>&x){return x*P[1]+P[2];},Data1);
        hist<> F2([&P](const value<>&x){return x*P[3]+P[4];},Data2);
	if(power==2){
	    F1+=[&P](const value<>&x)->value<>{return x*x*P[5];};
	    F2+=[&P](const value<>&x)->value<>{return x*x*P[6];};
	}
        const auto background1=toLine(F1);
        const auto background2=toLine(F2);
        const auto fit1=toLine(BW)*branching_ratio1.val()*P[0]+background1;
        const auto fit2=toLine(BW)*branching_ratio2.val()*P[0]+background2;
	const auto BINDING=static_cast<stringstream &>(stringstream()<< setprecision(4)<< B ).str();
	const auto GAMMA=static_cast<stringstream &>(stringstream()<< setprecision(4)<< G ).str();
        const string msg="B="+BINDING+" MeV; Г="+GAMMA+" MeV";
	string qua=(power==2)?"Quad":"";
	if((abs(la)>0.1)||(abs(lb)>0.1))qua=qua+"-"+to_string(int(la*10))+"_"+to_string(int(lb*10));
        Plot("UpperLimit-"+to_string(int(B*10))+"-"+to_string(int(G*10))+"Fit1"+qua,5)
            .Hist(Data1,"data").Line(fit1,"bg+signal","","lc rgb \"black\"").Line(background1,"bg","","lc rgb \"green\"")
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [10:40]"
                << "set title 'pd->^3He2γ "+msg+"'"<<"set key right top";
        Plot("UpperLimit-"+to_string(int(B*10))+"-"+to_string(int(G*10))+"Fit2"+qua,5)
            .Hist(Data2,"data").Line(fit2,"bg+signal","","lc rgb \"black\"").Line(background2,"bg","","lc rgb \"green\"")
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [15:50]"
                << "set title 'pd->^3He6γ "+msg+"'"<<"set key right top";
    }
    return FIT;
}
const double K=1.64485;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "bound-weighted-average");
    const auto binding=BinsByStep(-65.0,2.5,-2.5),gamma=BinsByStep(0.0,2.5,40.0);
    BiSortedPoints<value<>,value<>,Fitter> fit_results(binding,gamma,[](const ParamSet&){return 0.0;});
    cout<<"Fitting"<<endl;
    fit_results.FullCycleVar([](const value<>&B,const value<>&G,Fitter&res){
	res=BWfit("_",1,0,0,B.val(),G.val(),false);
    });
    BiSortedPoints<value<>,value<>,value<>> chisquare(binding,gamma),upper(binding,gamma);
    cout<<"converting"<<endl;
    Plot chisq_plot("UpperLimit-chisq-1d",4);
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
	    upper    .Bin(b,g)=A.val()+A.uncertainty()*K;
	}
	const auto WDTH=static_cast<stringstream &>(stringstream()<< setprecision(4)<< gamma[g].val()).str();
	if(g%4==3)chisq_plot.Line(chi_sq_curve,"Width="+WDTH+" MeV");
    }
    cout<<"plotting"<<endl;
    PlotHist2d(sp2,"UpperLimit-chisq").Distr(chisquare)<<"set colorbox"
	<<"set xlabel 'Peak position, MeV'"<<"set ylabel 'Width, MeV'"
	<<"set title 'χ^2/d'";
    PlotHist2d(sp2,"UpperLimit-upperlimit").Distr(upper)<<"set colorbox"
	<<"set xlabel 'Peak position, MeV'"<<"set ylabel 'Width, MeV'"
	<<"set title 'Upper limit, nb'";
    for(const auto&G:gamma){
        ext_hist<3> A_hist;
        SortedPoints<> upper_limit,lower_limit,systematics1,systematics2;
        for(const auto&B:binding){
	    RawSystematicError2 calc(params,[&B,&G](const string&suffix){
		function<Uncertainties<3>(const double&,const double&)> func=
			[&B,&G,&suffix](const double&power,const double&a)
		{
		    function<Uncertainties<3>(const double&)> F=[&B,&G,&suffix,&power,&a](const double&b){
			const auto fit=BWfit(suffix,size_t(power),a,b,B.val(),G.val(),
			    (suffix=="_")&&((size_t(power)==2)||(size_t(power)==1))&&
			    (G.Contains(9)||G.Contains(19)||G.Contains(29)||G.Contains(39))
			);
			const auto&A=fit.ParametersWithUncertainties()[0];
	        	return uncertainties(A.val(),A.uncertainty(),0.0,0.0);
		    };
		    return (suffix=="_")&&(power==1)&&(a==0)?SystematicError2<upper_limit_right>(F)(true):F(0);
		};
		return (suffix=="_")?SystematicError2<upper_limit_fit_power,upper_limit_left>(func)(true):func(1,0);
	    });
	    const auto A=calc(true);
	    A_hist<<make_point(B.val(),A);
	    upper_limit<<make_point(B.val(),A.val()+A.uncertainty<1>()*K);
	    lower_limit<<make_point(B.val(),(A.val()-A.uncertainty<1>()*K>0)?A.val()-A.uncertainty<1>()*K:0.0);
	    systematics1<<make_point(B.val(),(A.val()-A.uncertainty<2>()>0)?A.val()-A.uncertainty<2>():0.0);
	    systematics2<<make_point(B.val(),A.val()+A.uncertainty<3>());
        }
	const auto WDTH=static_cast<stringstream &>(stringstream()<< setprecision(4)<< G.val()).str();
	Plot("UpperLimit-Gconst"+to_string(int(G.val()*10)))
	    .Line(upper_limit,"Upper limit","","lc rgb \"blue\"")
	    .Line(systematics1,"Systematics","","lc rgb \"green\"")
	    .Line(systematics2,"Systematics","","lc rgb \"green\"")
	    <<"set title 'Width = "+WDTH+" MeV'"<<"set key on"
	    <<"set xlabel 'Peak position, MeV'"<<"set yrange [0:20]"
	    <<"set ylabel 'σ, nb'";
	Plot("UpperLimit-Gconst"+to_string(int(G.val()*10))+"-Parameter")
	    .Hist(take_uncertainty_component<1>(A_hist),"Fit resut")
	    .Line(upper_limit,"Stat. CL=90%","","lc rgb \"blue\"")
	    .Line(lower_limit,"","","lc rgb \"blue\"")
	    .Line(systematics1,"Systematics","","lc rgb \"green\"")
	    .Line(systematics2,"","","lc rgb \"green\"")
	    <<"set title 'Width = "+WDTH+" MeV'"<<"set key on"
	    <<"set xlabel 'Peak position, MeV'"<<"set yrange [-1:20]"
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
	Plot("UpperLimit-LinearFitChisq")
	    .Hist(chisq_hist,"","LinearFitChisq")
	    <<"set ylabel 'χ^2/d, n.d.'"
	    <<"set xlabel 'Parameter index'";
    }
}
