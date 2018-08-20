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
#include <Genetic/equation.h>
#include <Genetic/uncertainties.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Parameters/parameters.h>
#include <Parameters/systematic.h>
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

Chain<value_numeric_distr<>> BWfit(const string&suffix,const double&B,const double&G,bool plot=false){
    const auto Qbins=BinsByStep(-70.,2.5,2.5);
    const auto&lum=LuminosityGetter::Instance().get(suffix);
    cout<<"Histogram for "<<suffix<<": "<<"data1"<<endl;
    static Cache<string,ext_hist<2>> DATA1,DATA2,MBG1,MBG2;
    const auto&data1=DATA1(suffix,[&suffix,&lum,&Qbins](){
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
    cout<<"Histogram for "<<suffix<<": "<<"bg1"<<endl;
    const auto bg1=lum.XRange(-70,2.5)*MBG1(suffix,[&suffix,&Qbins](){
        return ext_hist<2>([&suffix](const value<>&Q){
            cout<<"Preparing histogram for "<<suffix<<": "<<"bg1["<<Q.min()<<":"<<Q.max()<<"]MeV"<<endl;
            const size_t bin_num=trunc((Q.val()+70.)/2.5);
            const auto MC_TIM=Hist(MC, "He3pi0pi0", histpath_central_reconstr1, "TIM6-Bin-"+to_string(bin_num),suffix).XRange(-0.3,0.3);
            static Cache<string,hist<>> NORM;
            const auto &N = NORM(suffix,[&suffix](){return Hist(MC, "He3pi0pi0", histpath_central_reconstr1, "0-Reference",suffix);})[bin_num].Y();
            return extend_value<2,2>(std_error(MC_TIM.TotalSum().val())/N);
        },Qbins);
    });
    cout<<"Histogram for "<<suffix<<": "<<"data2"<<endl;
    const auto&data2=DATA2(suffix,[&suffix,&lum,&Qbins](){
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
    cout<<"Histogram for "<<suffix<<": "<<"bg2"<<endl;
    const auto bg2=lum.XRange(-70,2.5)*MBG2(suffix,[&suffix,&Qbins](){
        return ext_hist<2>([&suffix](const value<>&Q){
            const size_t bin_num=trunc((Q.val()+70.)/2.5);
            cout<<"Preparing histogram for "<<suffix<<": "<<"bg2["<<Q.min()<<":"<<Q.max()<<"]MeV"<<endl;
            const auto MC_TIM=Hist(MC, "He3pi0pi0pi0", histpath_central_reconstr2, "TIM7-Bin-"+to_string(bin_num),suffix).XRange(-0.3,0.3);
            static Cache<string,hist<>> NORM;
            const auto &N = NORM(suffix,[&suffix](){return Hist(MC, "He3pi0pi0pi0", histpath_central_reconstr2, "0-Reference",suffix);})[bin_num].Y();
            return extend_value<2,2>(std_error(MC_TIM.TotalSum().val())/N);
        },Qbins);
    });
    const auto Data1=take_uncertainty_component<1>(data1.XRange(-70,5)),
               Data2=take_uncertainty_component<1>(data2.XRange(-70,5)),
               BG1=take_uncertainty_component<1>(bg1.XRange(-70,5)),
               BG2=take_uncertainty_component<1>(bg2.XRange(-70,5));
    cout<<"Fitting for "<<suffix<<": B="<<B<<"; G="<<G<<endl;
    const hist<> BW([&B,&G](const value<>&Q){
        //ToDo: calculate Breight Wigner using kinetic energy value
        return BreitWigner(Q.val(),B,G);
    },Data1);
    EquationSolver<DifferentialMutations<>,UncertaintiesEstimation> FIT{{
        .left=[&BW,&Data1,&BG1,&Data2,&BG2](const ParamSet&P){
            const auto F1=BW*branching_ratio1*P[0]+BG1*P[1];
            const auto F2=BW*branching_ratio2*P[0]+BG2*P[2];
            double res=0;
            for(size_t i=0;i<F1.size();i++){
                res+=F1[i].Y().NumCompare(Data1[i].Y());
                res+=F2[i].Y().NumCompare(Data2[i].Y());
            }
            return res;
        },
        .right=[](const ParamSet&){return 0;}
    }};
    FIT.SetFilter([](const ParamSet&P){return (P[0]>0)&&(P[1]>0)&&(P[2]>0);});
    FIT.SetUncertaintyCalcDeltas({0.01,0.01,0.01});
    FIT.Init(200,make_shared<InitialDistributions>()
        <<make_shared<DistribUniform>(0,30)
        <<make_shared<DistribUniform>(0,30000)
        <<make_shared<DistribUniform>(0,300)
    );
    while(!FIT.AbsoluteOptimalityExitCondition(0.000001))FIT.Iterate();
    cout << "Fit result: " << FIT.iteration_count() << " iterations; "
        << FIT.Optimality() << "<chi^2<"
        << FIT.Optimality(FIT.PopulationSize() - 1)
        << endl;
    cout << "xi^2/d = "<< FIT.Optimality()/(Data1.size()+Data2.size()-FIT.ParamCount())<<endl;
    const auto&P=FIT.Parameters();
    if(plot){
        const auto background1=toLine(BG1)*P[1];
        const auto background2=toLine(BG2)*P[2];
        const auto fit1=toLine((BW*branching_ratio1))*P[0]+background1;
        const auto fit2=toLine((BW*branching_ratio2))*P[0]+background2;
        const string msg="B="+to_string(B)+"; Г="+to_string(G);
        Plot("UpperLimitTotalFit1",5)
            .Hist(Data1).Line(fit1).Line(background1)
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:]"
                << "set title 'pd->3He+2g "+msg+"'"<<"set key right top";
        Plot("UpperLimitTotalFit2",5)
            .Hist(Data2).Line(fit2).Line(background2)
                <<"set key left top">>"set key right top"
                << "set xlabel 'Q, MeV'" << "set key on"<<"set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:]"
                << "set title 'pd->3He+6g "+msg+"'"<<"set key right top";
    }
    return FIT.ParametersWithUncertainties();
}
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "bound-weighted-average");
    const list<size_t> params{pbeam_corr,he3_cut_h,he3_theta_cut,
        gamma_E_thr,time_dt,time_t1,time_t2,eta_theta_thr,he3mm_cut,
        gamma_mm_lo,gamma_mm_hi,gamma_im_lo,gamma_im_hi,
        gamma_mm_lo,gamma_mm_hi,gamma_im_lo6,gamma_im_hi6,three_pi0
    };
    Plot plot_upper("UpperLimit",5);
    plot_upper
            <<"set key left top">>"set key right top"
            << "set xlabel 'Binding energy, MeV'" << "set key on"<<"set xrange [-70:10]"
            << "set ylabel 'peak height (fit), nb'" << "set yrange [-10:40]"
            << "set title 'Upper limit'"<<"set key right top";
    for(double G=5;G<=6;G+=10){
        ext_hist<2> UpperLimit;
        for(double B=-28.75;B<0;B+=2.5){
            UpperLimit<<make_point(B,RawSystematicError(params,[&G,&B](const string&suffix){
                return extend_value<1,2>(value<>(BWfit(suffix,B,G,(suffix=="_")&&(B>-10.)&&(B<-7.5)&&(G==5))[0]));
            })());
        }
        plot_upper.Hist_2bars<1,2>(UpperLimit,"Г="+to_string(G));
    }
}
