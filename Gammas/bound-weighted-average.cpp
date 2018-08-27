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
const value<> branching_ratio1{0.393,0.003};
const value<> branching_ratio2{0.322,0.003};
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
const auto Qbins=BinsByStep(-70.,2.5,0.0);
Fitter BWfit(const string&suffix,const double&B,const double&G,bool plot=false){
    const auto&lum=LuminosityGetter::Instance().get(suffix);
    cout<<"Histogram for "<<suffix<<": "<<"data1"<<endl;
    static Cache<string,ext_hist<2>> DATA1,DATA2;
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
    const hist<> BW([&B,&G](const value<>&Q){
	const auto b=Q2P(B),q=Q2P(Q.val());
	const auto g=G*0.001;
        return BreitWigner(q,b,g)/BreitWigner(b,b,g);
    },Data1);
    Fitter FIT([BW,Data1,Data2](const ParamSet&P){
        const auto F1=BW*branching_ratio1*P[0]+[&P](const value<>&x){return x*P[1]+P[2];};
        const auto F2=BW*branching_ratio2*P[0]+[&P](const value<>&x){return x*P[3]+P[4];};
        double res=0;
        for(size_t i=0;i<F1.size();i++){
            res+=F1[i].Y().NumCompare(Data1[i].Y());
            res+=F2[i].Y().NumCompare(Data2[i].Y());
        }
        return res;
    });
    FIT.SetFilter([](const ParamSet&P){return (P[0]>0);});
    FIT.Init(100,make_shared<InitialDistributions>()
        <<make_shared<DistribUniform>(0,10)
        <<make_shared<DistribUniform>(0,100)
        <<make_shared<DistribUniform>(0,10000)
        <<make_shared<DistribUniform>(0,100)
        <<make_shared<DistribUniform>(0,10000)
    );
    while(!FIT.AbsoluteOptimalityExitCondition(0.0000001))FIT.Iterate();
    cout << "Fit result: " << FIT.iteration_count() << " iterations; "
        << FIT.Optimality() << "<chi^2<"
        << FIT.Optimality(FIT.PopulationSize() - 1)
        << endl;
    FIT.SetUncertaintyCalcDeltas({0.1,0.1,0.01,0.1,0.01});
    const auto&P=FIT.ParametersWithUncertainties();
    if(plot){
        const hist<> background1([&P](const value<>&x){return x*P[1]+P[2];},BW);
        const hist<> background2([&P](const value<>&x){return x*P[3]+P[4];},BW);
        const auto fit1=BW*branching_ratio1*P[0]+background1;
        const auto fit2=BW*branching_ratio2*P[0]+background2;
        const string msg="B="+to_string(B)+"; Г="+to_string(G);
        Plot("UpperLimit-"+to_string(int(B*10))+"-"+to_string(int(G*10))+"Fit1",5)
            .Hist(Data1,"data").Line(toLine(fit1)).Line(toLine(background1))
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:]"
                << "set title 'pd->3He+2g "+msg+"'"<<"set key right top";
        Plot("UpperLimit-"+to_string(int(B*10))+"-"+to_string(int(G*10))+"Fit2",5)
            .Hist(Data2,"data").Line(toLine(fit2)).Line(toLine(background2))
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:]"
                << "set title 'pd->3He+6g "+msg+"'"<<"set key right top";
    }
    return FIT;
}
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "bound-weighted-average");
    const list<size_t> params{pbeam_corr,he3_cut_h,he3_theta_cut,
        gamma_E_thr,time_dt,time_t1,time_t2,eta_theta_thr,he3mm_cut,
        gamma_mm_lo,gamma_mm_hi,gamma_im_lo,gamma_im_hi,
        gamma_mm_lo,gamma_mm_hi,gamma_im_lo6,gamma_im_hi6,three_pi0
    };
    const auto binding=BinsByStep(-60.0,2.5,0.0),gamma=BinsByStep(0.0,2.5,60.0);
    BiSortedPoints<value<>,value<>,Fitter> fit_results(binding,gamma,[](const ParamSet&){return 0.0;});
    cout<<"Fitting"<<endl;
    fit_results.FullCycleVar([](const value<>&B,const value<>&G,Fitter&res){
	res=BWfit("_",B.val(),G.val(),G.Contains(19));
    });
    BiSortedPoints<value<>,value<>,value<>> chisquare(binding,gamma),upper(binding,gamma);
    cout<<"converting"<<endl;
    for(size_t b=0;b<binding.size();b++)for(size_t g=0;g<gamma.size();g++){
	const auto&src=fit_results[b][g];
	chisquare.Bin(b,g)=src.Optimality()/(Qbins.size()+Qbins.size()-src.ParamCount());
	upper    .Bin(b,g)=src.ParametersWithUncertainties()[0].max();
    }
    cout<<"plotting"<<endl;
    PlotHist2d(sp2,"UpperLimit-chisq").Distr(chisquare)<<"set colorbox"
	<<"set xlabel 'Peak position, MeV'"<<"set ylabel 'Width, MeV'"
	<<"set title 'Chi square'";
    PlotHist2d(sp2,"UpperLimit-upperlimit").Distr(upper)<<"set colorbox"
	<<"set xlabel 'Peak position, MeV'"<<"set ylabel 'Width, MeV'"
	<<"set title 'Upper limit, nb'";
    ext_hist<2> UpperLimitMore;
    for(const auto&B:binding){
	const double G=18.75;
	UpperLimitMore<<make_point(B.val(),RawSystematicError(params,[&B,&G](const string&suffix){
	    const auto fit=BWfit(suffix,B.val(),G,false);
            return uncertainties(fit.ParametersWithUncertainties()[0].max(),0.0,0.0);
        })());	
    }
    Plot("UpperLimit-Gconst")
	.Hist(wrap_hist(UpperLimitMore))
	<<"set title 'G = 18.75 MeV'"
	<<"set xlabel 'Peak position, MeV'"
	<<"set ylabel 'Upper limit, nb'";
}