// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <math_h/sigma3.h>
#include <Genetic/fit.h>
#include <Genetic/paramfunc.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Kinematics/particles.h>
using namespace std;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
double Q2E(const double&Q){//The Breight-Wigner fotmula takes not Q but projectile energy
    const auto CM=binaryDecay(Particle::he3().mass()+Particle::eta().mass()+Q*0.001,Particle::p().mass(),Particle::d().mass());
    const auto p_lab=CM.first.Transform(CM.second.Beta());
    return p_lab.Ekin()*1000.;
}
typedef Mul<Par<0>,Func3<BreitWigner,Func<Q2E,Arg<0>>,Func<Q2E,Par<1>>,Div<Par<2>,Const<2>>>> FG;
typedef PolynomFunc<Arg<0>,3,1> BG;
typedef PolynomFunc<Arg<0>,0,1> BG2;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma");
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound1-2g","bound2-2g","bound3-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
    const vector<string> suffix={"-m60","-m40","-m20","-0"};
    const ext_hist<2> luminosity = Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc");
    const ext_hist<2> luminosity_he = Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYf");
    const auto true_he3eta = luminosity_he
        *extend_hist<2,2>(hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed")))
            .XRange(luminosity_he.left().X().min(),luminosity_he.right().X().max())
        /trigger_he3_forward.scaling;
    const double branching_ratio=0.39;
    hist<> CS,POS,WIDTH;
    SortedPoints<> CHISQ,CHISQ_W;
    for(size_t a_t=0;a_t<suffix.size();a_t++){
        cout<<suffix[a_t]<< " fitting"<<endl;
        const double cutpos=-0.06+0.02*a_t;
        const hist<> acc_bound=Plotter::Instance().GetPoints<value<>>("He3gg-acceptance"+suffix[a_t]+"-2");
        const hist<> acc_he3eta=Plotter::Instance().GetPoints<value<>>("He3gg-acceptance"+suffix[a_t]+"-3");
        cout<<suffix[a_t]<< " fitting"<<endl;
        const hist<> ev=Plotter::Instance().GetPoints<value<>>("He3gg-data"+suffix[a_t]);
        const auto known_events = (true_he3eta*branching_ratio)
            *extend_hist<2,2>(acc_he3eta).XRange(true_he3eta.left().X().min(),true_he3eta.right().X().max());
        cout<<suffix[a_t]<< " fitting"<<endl;
        const auto data_shape=(
            extend_hist<1,2>(ev)*trigger_he3_forward.scaling/(extend_hist<2,2>(acc_bound)*luminosity)
        ).XRange(-60,2.5);
        FitFunction2<DifferentialMutations<>,Add<FG,BG>> fit(wrap_hist(data_shape).removeXerorbars());
        auto init=make_shared<InitialDistributions>()
                    <<make_shared<DistribGauss>(50,50)
                    <<make_shared<DistribGauss>(-15,5)
                    <<make_shared<DistribGauss>(10,10)
                    <<make_shared<DistribGauss>(100,100);
        while(init->Count()<BG::ParamCount)init<<make_shared<DistribGauss>(0,1);
        fit.SetFilter([](const ParamSet&P){
            return (P[0]>0)&&(P[1]<0)&&(P[1]>-30)&&(P[2]>9)&&(P[2]<40);
        });
        fit.Init(500,init);
        while(!fit.AbsoluteOptimalityExitCondition(0.0000001))fit.Iterate();
        fit.SetUncertaintyCalcDeltas(parEq(BG::ParamCount,0.01));
        const auto&P=fit.ParametersWithUncertainties();
        const auto chain=ChainWithStep(-70.,0.01,2.5);
        const SortedPoints<>
        fg([&fit](double x){return fit({x});},chain),
        bg([&fit](double x){return BG()({x},fit.Parameters());},chain);
        Plot("He3gg-events-norm"+suffix[a_t]+"-bound")
            .Hist_2bars<1,2>(data_shape,"Data")
            .Line(bg).Line(fg,"fit with peak")
                << "set xlabel 'Q, MeV'" << "set key on"
                << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]";
        //cross section is not peak area but it's height
        const value_numeric_const<> cs_coef=LinearInterpolation<>(fg-bg)(P[1].val())/P[0].val();
        CS<<make_point(cutpos*1000,P[0]*cs_coef);
        POS<<make_point(cutpos*1000,P[1]);
        WIDTH<<make_point(cutpos*1000,P[2]);
        CHISQ<<make_point(cutpos*1000,fit.Optimality()/(fit.Points().size()-fit.ParamCount()));
        //pure background fit
        FitFunction<DifferentialMutations<>,BG2> fit_w(wrap_hist(data_shape).removeXerorbars());
        auto init_w=make_shared<InitialDistributions>()
                <<make_shared<DistribGauss>(100,100);
        while(init_w->Count()<BG2::ParamCount)init_w<<make_shared<DistribGauss>(0,1);
        fit_w.Init(300,init_w);
        while(!fit_w.AbsoluteOptimalityExitCondition(0.0000001))fit_w.Iterate();
        CHISQ_W<<make_point(cutpos*1000,fit_w.Optimality()/(fit_w.Points().size()-fit_w.ParamCount()+3));
        const SortedPoints<> bg2([&fit_w](double x){return fit_w({x});},chain);
        Plot("He3gg-events-norm2"+suffix[a_t]+"-bound")
            .Hist_2bars<1,2>(data_shape,"Data")
            .Line(bg2,"fit without peak")
                << "set xlabel 'Q, MeV'" << "set key on"
                << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]";
    }
    cout<<"Final plots"<<endl;
    Plot("He3gg-cross-section").Hist(CS,"data")
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-70:30]"<<"set key on"
            << "set ylabel 'Cross section, nb'" << "set yrange [0:30]"
            << "set title 'Cross section estimation'";
    Plot("He3gg-cross-section-2").Hist(CS/branching_ratio,"divided by branching ratio")
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-70:30]"<<"set key on"
            << "set ylabel 'Cross section, nb'" << "set yrange [0:]"
            << "set title 'Cross section estimation'";
    Plot("He3gg-pos").Hist(POS)
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-70:30]"
            << "set ylabel 'Position, MeV'" << "set yrange [-40:0]"
            << "set title 'Peak position'";
    Plot("He3gg-width").Hist(WIDTH)
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-70:30]"
            << "set ylabel 'Gamma, MeV'" << "set yrange [0:]"
            << "set title 'Peak width (Gamma)'";
    Plot("He3gg-cross-section-chisq").Line(CHISQ)
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-70:30]"<<"set key on"
            << "set ylabel 'chi square, n.d.'" << "set yrange [0:2]"
            << "set title 'Chi square'";
    Plot("He3gg-cross-section-chisq2").Line(CHISQ,"with peak").Line(CHISQ_W,"without peak")
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-70:30]"<<"set key on"
            << "set ylabel 'chi square, n.d.'" << "set yrange [0:2]"
            << "set title 'Chi square'";
    cout<<"END"<<endl;
}

