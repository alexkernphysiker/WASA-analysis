// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Genetic/fit.h>
#include <Genetic/paramfunc.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
typedef Mul<Par<0>,Func3<BreitWigner,Arg<0>,Par<1>,Par<2>>> FG;
typedef PolynomFunc<Arg<0>,3,2> BG;
typedef PolynomFunc<Arg<0>,0,2> BG2;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma");
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound3-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("C");
    const hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";

    const vector<string> suffix={"-m40","-m20","-0","-20","-40"};
    const hist<> luminosity = Plotter::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> luminosity_he = Plotter::Instance().GetPoints<value<>>("LUMINOSITYf");
    const hist<> true_he3eta = luminosity_he
        *hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed"))
            .XRange(luminosity_he.left().X().min(),luminosity_he.right().X().max())
        /trigger_he3_forward.scaling;
    const double branching_ratio=0.39;
    hist<> CS,POS,WIDTH;
    SortedPoints<> CHISQ,CHISQ_W;
    for(size_t a_t=0;a_t<suffix.size();a_t++){
        cout<<suffix[a_t]<< " fitting"<<endl;
        const auto TIM=Hist(DATA, "C", histpath_central_reconstr,"TIM4-AllBins");
        const double high=TIM.TransponateAndSort().right().X().max();
        const double cutpos=-0.04+0.02*a_t;
        Plot("He3gg-above-tim-data-allbins"+suffix[a_t]).Hist(TIM)
        .Line({make_point(cutpos,0.),make_point(cutpos,high)})
            << "set key on" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";

        const hist<> acc_bound=Plotter::Instance().GetPoints<value<>>("He3gg-acceptance"+suffix[a_t]+"-0");
        const hist<> acc_he3eta=Plotter::Instance().GetPoints<value<>>("He3gg-acceptance"+suffix[a_t]+"-1");
        cout<<suffix[a_t]<< " fitting"<<endl;
        const hist<> ev=Plotter::Instance().GetPoints<value<>>("He3gg-data"+suffix[a_t]);
        const hist<> known_events = (true_he3eta*branching_ratio)
            *acc_he3eta.XRange(true_he3eta.left().X().min(),true_he3eta.right().X().max());
        Plot("He3gg-events-final"+suffix[a_t]+"-bound")
            .Hist(ev.XRange(-50,0),"data, below threshold")
            .Hist(ev.XRange(12.5,30)-known_events,"data-3Heeta, upper threshold")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-45:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
        cout<<suffix[a_t]<< " fitting"<<endl;
        const auto data_shape=(
            (ev*trigger_he3_forward.scaling)/(acc_bound*luminosity)
        ).XRange(-50,0);
        FitFunction2<DifferentialMutations<>,Add<FG,BG>> fit(data_shape.removeXerorbars());
        auto init=make_shared<InitialDistributions>()
                    <<make_shared<DistribGauss>(50,50)
                    <<make_shared<DistribGauss>(-15,5)
                    <<make_shared<DistribGauss>(10,10)
                    <<make_shared<DistribGauss>(100,100);
        while(init->Count()<BG::ParamCount)init<<make_shared<DistribGauss>(0,1);
        fit.SetFilter([](const ParamSet&P){
            return (P[0]>0)&&(P[1]<0)&&(P[1]>-30)&&(P[2]>2.5)&&(P[2]<20);
        });
        fit.Init(300,init);
        while(!fit.AbsoluteOptimalityExitCondition(0.0000001))fit.Iterate();
        fit.SetUncertaintyCalcDeltas({0.1,0.01,0.01,0.1});
        const auto&P=fit.ParametersWithUncertainties();
        const auto chain=ChainWithStep(-50.,0.01,0.);
        const SortedPoints<>
        fg([&fit](double x){return fit({x});},chain),
        bg([&fit](double x){return BG()({x},fit.Parameters());},chain);
        Plot("He3gg-events-norm"+suffix[a_t]+"-bound")
            .Hist(data_shape,"Data")
            .Line(bg).Line(fg,"fit with peak")
                << "set xlabel 'Q, MeV'" << "set key on"
                << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
        //cross section is not peak area but it's height
        const value_numeric_const<> cs_coef=LinearInterpolation<>(fg-bg)(P[1].val())/P[0].val();
        CS<<make_point(cutpos*1000,P[0]*cs_coef);
        POS<<make_point(cutpos*1000,P[1]);
        WIDTH<<make_point(cutpos*1000,P[2]);
        CHISQ<<make_point(cutpos*1000,fit.Optimality()/(fit.Points().size()-fit.ParamCount()));
        //pure background fit
        FitFunction<DifferentialMutations<>,BG2> fit_w(data_shape.removeXerorbars());
        auto init_w=make_shared<InitialDistributions>()
                <<make_shared<DistribGauss>(100,100);
        while(init_w->Count()<BG2::ParamCount)init_w<<make_shared<DistribGauss>(0,1);
        fit_w.Init(300,init_w);
        while(!fit_w.AbsoluteOptimalityExitCondition(0.0000001))fit_w.Iterate();
        CHISQ_W<<make_point(cutpos*1000,fit_w.Optimality()/(fit_w.Points().size()-fit_w.ParamCount()+3));
        const SortedPoints<> bg2([&fit_w](double x){return fit_w({x});},chain);
        Plot("He3gg-events-norm2"+suffix[a_t]+"-bound")
            .Hist(data_shape,"Data")
            .Line(bg2,"fit without peak")
                << "set xlabel 'Q, MeV'" << "set key on"
                << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
    }
    cout<<"Final plots"<<endl;
    Plot("He3gg-cross-section").Hist(CS,"data")
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-50:50]"<<"set key on"
            << "set ylabel 'Cross section, nb'" << "set yrange [0:]"
            << "set title 'Cross section "+runmsg+"'";
    Plot("He3gg-cross-section-2").Hist(CS/branching_ratio,"divided by branching ratio")
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-50:50]"<<"set key on"
            << "set ylabel 'Cross section, nb'" << "set yrange [0:]"
            << "set title 'Cross section "+runmsg+"'";
    Plot("He3gg-pos").Hist(POS)
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-50:50]"
            << "set ylabel 'Position, MeV'" << "set yrange [-20:0]"
            << "set title 'Peak position "+runmsg+"'";
    Plot("He3gg-width").Hist(WIDTH)
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-50:50]"
            << "set ylabel 'sigma, MeV'" << "set yrange [0:20]"
            << "set title 'Peak width (sigma) "+runmsg+"'";
    Plot("He3gg-cross-section-chisq").Line(CHISQ)
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-50:50]"<<"set key on"
            << "set ylabel 'chi square, n.d.'" << "set yrange [0:2]"
            << "set title 'Chi square "+runmsg+"'";
    Plot("He3gg-cross-section-chisq2").Line(CHISQ,"with peak").Line(CHISQ_W,"without peak")
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-50:50]"<<"set key on"
            << "set ylabel 'chi square, n.d.'" << "set yrange [0:]"
            << "set title 'Chi square "+runmsg+"'";
    cout<<"END"<<endl;
}

