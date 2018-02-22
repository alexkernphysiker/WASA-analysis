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
typedef PolynomFunc<Arg<0>,3,3> BG;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma");
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound3-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("C");
    const hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";

    Plot gE("He3gg-gamma-energy-theory"), gEc("He3gg-gamma-energy-theory-cut"),
         gEd("He3gg-gamma-energy-data");
    gE << "set key on" << "set title 'Gamma energy. MC'"
       << "set yrange [0:]" << "set xlabel 'gamma energy, GeV'";
    gEc << "set key on" << "set title 'Gamma energy. MC. Cut'"
        << "set yrange [0:]" << "set xlabel 'gamma energy, GeV'";
    gEd << "set key on" << "set title 'Gamma energy. Data'"
        << "set yrange [0:]" << "set xlabel 'gamma energy, GeV'";
    for (size_t i = 0; i < reaction.size(); i++) {
        const auto &r = reaction[i];
        gE.Hist(
            Hist(MC, r, histpath_central_reconstr, "GammaEnergy")
            / Hist(MC, r, histpath_central_reconstr, "0-Reference").TotalSum().val()
            , r);
        gEc.Hist(
            Hist(MC, r, histpath_central_reconstr, "GammaEnergyCut")
            / Hist(MC, r, histpath_central_reconstr, "0-Reference").TotalSum().val()
            , r);
    }
    gEd.Hist(Hist(DATA, "C", histpath_central_reconstr, "GammaEnergy")).Hist(Hist(DATA, "C", histpath_central_reconstr, "GammaEnergyCut"));
    Plot("He3gg-gamma-count").Hist(Hist(DATA, "C", histpath_central_reconstr, "GammaCount"))<<"set log y">>"unset log y";

    vector<hist<>> ev_am;
    vector<vector<hist<>>> acceptance;
    hist<> he3acc;

    const vector<string> suffix={"-0","-20","-40"};
    vector<hist<>> ev_am1,ev_am2;
    vector<vector<hist<>>> acc;
    for(size_t i=0;i<suffix.size();i++){
        acc.push_back({});
        for (size_t j = 0; j < reaction.size(); j++)
            acc[i].push_back(hist<>());
        ev_am1.push_back(hist<>());
        ev_am2.push_back(hist<>());
    }

    const size_t cuts_count=10;
    for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
        acceptance.push_back(vector<hist<>>());
        for (size_t i = 0; i < reaction.size(); i++) {
            acceptance[cut_index].push_back(hist<>());
        }
        ev_am.push_back(hist<>());
    }
    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
                            << "Q in [" << setprecision(3)
                            << Q.min() << "; " << Q.max() << "] MeV").str();
        {
            hist<> he3reg = Hist(MC, reaction[0], histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num));
            he3acc<<make_point(Q,std_error(he3reg.TotalSum().val())/norm[bin_num].Y());
        }
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            Plot(
                Q.Contains(21) ? "He3gg-above-he3mm-mc-" + r : (
                    Q.Contains(-39) ? "He3gg-below-he3mm-mc-" + r : (
                        Q.Contains(-11) ? "He3gg-thr-he3mm-mc-" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+2gamma")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggmm-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-ggmm-mc" + r : (
                        Q.Contains(-11) ? "He3gg-thr-ggmm-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggim-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-ggim-mc" + r : (
                        Q.Contains(-11) ? "He3gg-thr-ggim-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-tim-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-tim-mc" + r : (
                        Q.Contains(-11) ? "He3gg-thr-tim-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-t-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-t-mc"+r : (
                        Q.Contains(-11) ? "He3gg-thr-t-mc"+r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("t4-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-dt-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-dt-mc"+r : (
                        Q.Contains(-11) ? "He3gg-thr-dt-mc"+r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("dt4-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";

        }{
            Plot(
                Q.Contains(21) ? "He3gg-above-he3mm-data" : (
                    Q.Contains(-39) ? "He3gg-below-he3mm-data" : (
                        Q.Contains(-11) ? "He3gg-thr-he3mm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+2gamma")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggmm-data" : (
                    Q.Contains(-39) ? "He3gg-below-ggmm-data" : (
                        Q.Contains(-11) ? "He3gg-thr-ggmm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggim-data" : (
                    Q.Contains(-39) ? "He3gg-below-ggim-data" : (
                        Q.Contains(-11) ? "He3gg-thr-ggim-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-tim-data" : (
                    Q.Contains(-39) ? "He3gg-below-tim-data" : (
                        Q.Contains(-11) ? "He3gg-thr-tim-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";

            Plot(
                Q.Contains(21) ? "He3gg-above-data-t" : (
                    Q.Contains(-39) ? "He3gg-below-data-t" : (
                        Q.Contains(-11) ? "He3gg-thr-data-t" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("t4-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-data-dt" : (
                    Q.Contains(-39) ? "He3gg-below-data-dt" : (
                        Q.Contains(-11) ? "He3gg-thr-data-dt" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("dt4-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
        }
        hist<> ACC,ACCbg,RATIO;
        for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
            const double TIM_DIFF_CUT=0.01*cut_index;
            value<> R=INFINITY;
            for (size_t i = 0; i < reaction.size(); i++) {
                const auto &r = reaction[i];
                hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
                const auto &N = Norm[bin_num].Y();
                if (N.Above(0)) {
                    const hist<> h = Hist(MC, r, histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)).XRange(TIM_DIFF_CUT,INFINITY);
                    const auto a=std_error(h.TotalSum().val()) / N;
                    acceptance[cut_index][i] << point<value<>>(Q,a);
                    if(i==0){
                        R=a;
                        ACC<<make_point(10.*cut_index,a);
                    }
                    if(i==2){
                        R/=a;
                        ACCbg<<make_point(10.*cut_index,a);
                    }
                } else {
                    acceptance[cut_index][i] << point<value<>>(Q, 0.0);
                    if(i==0)R=0.0;
                    if(i==2)R=INFINITY;
                }
            }
            if(isfinite(R.val()))RATIO<<make_point(0.001*cut_index,R);
            const hist<> data = Hist(DATA, "C", histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)).XRange(TIM_DIFF_CUT,INFINITY);
            ev_am[cut_index] << point<value<>>(Q, std_error(data.TotalSum().val()));
        }
        if(RATIO.size()>0)
        Plot(
            Q.Contains(21) ? "He3gg-above-acceptance-ratio" : (
                Q.Contains(-39) ? "He3gg-below-acceptance-ratio" : (
                    Q.Contains(-11) ? "He3gg-thr-acceptance-ratio" : ""
                )
            )
        )
        .Hist(RATIO)
                << "set key on" << "set title 'Acceptance ratio" + Qmsg + "'"
                << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'";
        if(ACC.size()>0)
        Plot(
            Q.Contains(21) ? "He3gg-above-acceptance" : (
                Q.Contains(-39) ? "He3gg-below-acceptance" : (
                    Q.Contains(-11) ? "He3gg-thr-acceptance" : ""
                )
            )
        )
        .Hist(ACC*100)
                << "set key on" << "set title 'Acceptance for foreground " + Qmsg + "'"
                << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
                << "set ylabel 'Acceptance, percents'"<<"set xrange [-5:100]";
        if(ACCbg.size()>0)
        Plot(
            Q.Contains(21) ? "He3gg-above-acceptance-bg" : (
                Q.Contains(-39) ? "He3gg-below-acceptance-bg" : (
                    Q.Contains(-11) ? "He3gg-thr-acceptance-bg" : ""
                )
            )
        )
        .Hist(ACCbg*100)
                << "set key on" << "set title 'Acceptance for the most intensive background " + Qmsg + "'"
                << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
                << "set ylabel 'Acceptance, percents'"<<"set xrange [-5:100]";

        for(size_t i = 0; i < reaction.size(); i++) for(size_t a_t=0;a_t<suffix.size();a_t++){
            const auto &r = reaction[i];
            const auto DT=Hist(MC, r, histpath_central_reconstr, "dt5"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(10);
            const auto T=Hist(MC, r, histpath_central_reconstr, "t5"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(10);
            hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) acc[a_t][i] << make_point(Q, std_error(DT.TotalSum().val())/N);
            else acc[a_t][i] << make_point(Q, 0.0);
            Plot(
                Q.Contains(21) ? "He3gg-above-dt-final"+suffix[a_t]+"-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-dt-final"+suffix[a_t]+"-mc"+r : (
                        Q.Contains(-11) ? "He3gg-thr-dt-final"+suffix[a_t]+"-mc"+r : ""
                    )
                )
            )
            .Hist(DT)
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-t-final"+suffix[a_t]+"-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-t-final"+suffix[a_t]+"-mc"+r : (
                        Q.Contains(-11) ? "He3gg-thr-t-final"+suffix[a_t]+"-mc"+r : ""
                    )
                )
            )
            .Hist(T)
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
        }
        for(size_t a_t=0;a_t<suffix.size();a_t++){
            const auto DT=Hist(DATA, "C", histpath_central_reconstr, "dt5"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(50);
            const auto T=Hist(DATA, "C", histpath_central_reconstr, "t5"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(50);

            const auto DTBG=WeightedAverage<>()<<DT[2].Y()<<DT[3].Y()<<DT[4].Y();
            ev_am1[a_t]<<make_point(Q,DT[0].Y()+DT[1].Y()-DTBG*2.0);
            const LinearInterpolation<value<>> TBG=Points<value<>>{T[3],T[7]};
            ev_am2[a_t]<<make_point(Q,(T[4].Y()-TBG(T[4].X()))+(T[5].Y()-TBG(T[5].X()))+(T[6].Y()-TBG(T[6].X())));

            const hist<> dtbgplot=Points<value<>>{{DT[0].X(),DTBG},{DT[1].X(),DTBG}};
            Plot(
                Q.Contains(21) ? "He3gg-above-data-dt-final"+suffix[a_t] : (
                    Q.Contains(-39) ? "He3gg-below-data-dt-final"+suffix[a_t] : (
                        Q.Contains(-11) ? "He3gg-thr-data-dt-final"+suffix[a_t] : ""
                    )
                )
            )
            .Hist(DT,"Data").Hist(dtbgplot,"background")
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<<"set key on";
            const hist<> tbgplot=Points<value<>>{{T[4].X(),TBG(T[4].X())},{T[5].X(),TBG(T[5].X())},{T[6].X(),TBG(T[6].X())}};
            Plot(
                Q.Contains(21) ? "He3gg-above-data-t-final"+suffix[a_t] : (
                    Q.Contains(-39) ? "He3gg-below-data-t-final"+suffix[a_t] : (
                        Q.Contains(-11) ? "He3gg-thr-data-t-final"+suffix[a_t] : ""
                    )
                )
            )
            .Hist(T,"Data").Hist(tbgplot,"background")
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'"<<"set key on";
        }

    }
    const hist<> luminosity = Plotter::Instance().GetPoints<value<>>("LUMINOSITYc");
    hist<> true_he3eta = Plotter::Instance().GetPoints<value<>>("LUMINOSITYf");
    true_he3eta*=Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed");
    true_he3eta/= double(trigger_he3_forward.scaling);
    while(true_he3eta.left().X().min()>luminosity.left().X().min())
        true_he3eta<<make_point(value<>(true_he3eta.left().X().min()-1.25,1.25),0);
    for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
        Plot accplot("He3gg-acceptance-"+to_string(cut_index));
        accplot << "set title 'Cut number "+to_string(cut_index)+"'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, percents'"
            << "set yrange [0.0001:500]" << "set xrange [-70:30]"
            << "set key on"<<"set log y">>"unset log y";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto acc = acceptance[cut_index][i].YRange(0.0000001, INFINITY);
            if (acc.size() > 0) {
                accplot.Hist(acc*100, reaction[i]);
            }
        }
        Plot("He3gg-events-"+to_string(cut_index))
        .Hist(ev_am[cut_index],"All data")
            << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title 'Cut number "+to_string(cut_index)+". "+runmsg+"'";
    }
    hist<> CS_1,CS_3,POS,WIDTH;
    SortedPoints<> CHISQ,CHISQ_W;
    for(size_t a_t=0;a_t<suffix.size();a_t++){
        Plot accplot("He3gg-acceptance-final"+suffix[a_t]);
        accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, percents'"
            << "set yrange [0.00001:500]" << "set xrange [-70:30]"
            << "set key on"<<"set log y">>"unset log y";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto ac = acc[a_t][i].YRange(0.0000001, INFINITY);
            if (ac.size() > 0) {
                accplot.Hist(ac*100, reaction[i]);
                Plot("He3gg-acceptance-final"+suffix[a_t]+"-"+reaction[i]).Hist(ac*100)
                    << "set title 'Acceptance "+reaction[i]+"'"
                    << "set xlabel 'Q, MeV'"
                    << "set ylabel 'Acceptance, percents'"
                    << "set xrange [-70:30]";
            }
        }
        const hist<> known_events = true_he3eta*acc[a_t][1]*0.4;
        Plot("He3gg-events-final"+suffix[a_t])
            .Hist(ev_am1[a_t],"data dt(gamma-gamma)")
            .Hist(ev_am2[a_t],"data dt(3He-gamma)")
            .Hist(known_events,"3He+eta")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
        const auto ev=hist_avr(ev_am1[a_t],ev_am2[a_t]);
        Plot("He3gg-events-final"+suffix[a_t]+"-bound")
            .Hist(ev - known_events)
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
        const auto data_shape=(
            (ev-known_events)*double(trigger_he3_forward.scaling)/(acc[a_t][0]*luminosity)
        ).XRange(-50,30);
        FitFunction2<DifferentialMutations<>,Add<FG,BG>> fit(data_shape.removeXerorbars());
        auto init=make_shared<InitialDistributions>()
                    <<make_shared<DistribGauss>(50,50)
                    <<make_shared<DistribGauss>(-15,5)
                    <<make_shared<DistribGauss>(10,10)
                    <<make_shared<DistribGauss>(100,100);
        while(init->Count()<BG::ParamCount)init<<make_shared<DistribGauss>(0,1);
        fit.SetFilter([](const ParamSet&P){return (P[2]>4)&&(P[2]<25)&&(P[0]>0)&&(P[1]<5)&&(P[1]>-30);});
        fit.Init(300,init);
        while(!fit.AbsoluteOptimalityExitCondition(0.0000001))fit.Iterate();
        fit.SetUncertaintyCalcDeltas({0.1,0.01,0.01,0.1});
        const auto&P=fit.ParametersWithUncertainties();
        const auto chain=ChainWithStep(-50.,0.001,30.);
        const SortedPoints<>
        fg([&fit](double x){return fit({x});},chain),
        bg([&fit](double x){return BG()({x},fit.Parameters());},chain);
        Plot("He3gg-events-norm"+suffix[a_t]+"-bound")
            .Hist(data_shape,"Data")
            .Line(bg).Line(fg)
                << "set xlabel 'Q, MeV'" << "set key on"
                << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
        //cross section is not peak area but it's height
        const value_numeric_const<> cs_coef=LinearInterpolation<>(fg-bg)(P[1].val())/P[0].val();
        const double cp=20.*a_t;
        CS_1<<make_point(cp,P[0]*cs_coef);
        CS_3<<make_point(cp,value<>(P[0]).make_wider(3)*cs_coef);
        //other parameters
        POS<<make_point(cp,P[1]);
        WIDTH<<make_point(cp,P[2]);
        CHISQ<<make_point(cp,fit.Optimality()/(fit.Points().size()-fit.ParamCount()));
        //pure background fit
        FitFunction<DifferentialMutations<>,BG> fit_w(data_shape.removeXerorbars());
        auto init_w=make_shared<InitialDistributions>()
                <<make_shared<FixParam>(0)
                <<make_shared<FixParam>(0)
                <<make_shared<FixParam>(0)
                <<make_shared<DistribGauss>(100,100);
        while(init_w->Count()<BG::ParamCount)init_w<<make_shared<DistribGauss>(0,1);
        fit_w.Init(300,init_w);
        while(!fit_w.AbsoluteOptimalityExitCondition(0.0000001))fit_w.Iterate();
        CHISQ_W<<make_point(cp,fit_w.Optimality()/(fit_w.Points().size()-fit_w.ParamCount()+3));
    }
    Plot("He3gg-tube-acc").Hist(he3acc*100.)
            << "set xlabel 'Q, MeV'" << "set xrange [-70:30]"
            << "set ylabel 'Acceptance, percents'" << "set yrange [0:]"
            << "set title 'How many helium ions from mesic nuclei decay would be detected'";
    Plot("He3gg-cross-section").Hist(CS_3,"3 sigma").Hist(CS_1,"1 sigma")
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-10:50]"<<"set key on"
            << "set ylabel 'Cross section, nb'" << "set yrange [0:]"
            << "set title 'Cross section "+runmsg+"'";
    Plot("He3gg-pos").Hist(POS)
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-10:50]"
            << "set ylabel 'Position, MeV'" << "set yrange [-20:0]"
            << "set title 'Peak position "+runmsg+"'";
    Plot("He3gg-width").Hist(WIDTH)
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-10:50]"
            << "set ylabel 'sigma, MeV'" << "set yrange [0:20]"
            << "set title 'Peak width (sigma) "+runmsg+"'";
    Plot("He3gg-cross-section-chisq").Line(CHISQ,"with peak").Line(CHISQ_W,"without peak")
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, MeV'"
            << "set xrange [-10:50]"<<"set key on"
            << "set ylabel 'chi square, n.d.'" << "set yrange [0:2]"
            << "set title 'Chi square "+runmsg+"'";
}
