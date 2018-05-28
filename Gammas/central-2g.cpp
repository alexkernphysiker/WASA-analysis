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
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma");
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound1-2g","bound2-2g","bound3-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("All");
    const hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    cout<<"First plots"<<endl;
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
    gEd.Hist(Hist(DATA, "All", histpath_central_reconstr, "GammaEnergy")).Hist(Hist(DATA, "All", histpath_central_reconstr, "GammaEnergyCut"));
    Plot("He3gg-gamma-count").Hist(Hist(DATA, "All", histpath_central_reconstr, "GammaCount"));

    vector<hist<>> ev_am;
    vector<vector<hist<>>> acceptance;
    const vector<string> suffix={"-m60","-m40","-m20","-0","-20","-40"};
    vector<vector<hist<>>> acc;
    for(size_t i=0;i<suffix.size();i++){
        acc.push_back({});
        for (size_t j = 0; j < reaction.size(); j++)
            acc[i].push_back(hist<>());
        ev_am.push_back(hist<>());
    }

    for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            cout<<"All-bins MC plots "<<r<<endl;
            Plot("He3gg-he3mm-mc-" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM0"))
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM1"), "cut")
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM2"), "2 gammas required")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [0.45:0.57]"
                    << "set xlabel '3He missing mass - Q, GeV'";

            PlotHist2d(sp2, "He3gg-alpha1-mc" + r)
                .Distr(Hist2d(MC,r,histpath_central_reconstr,"GGangle2").Scale(5,5))
                   << "set key on"
                   << "set xlabel 'sin alpha'"
                   << "set ylabel 'cos alpha'";
            PlotHist2d(sp2, "He3gg-alpha2-mc" + r)
                .Distr(Hist2d(MC,r,histpath_central_reconstr,"GGangle3").Scale(5,5))
                   << "set key on"
                   << "set xlabel 'sin alpha'"
                   << "set ylabel 'cos alpha'";

            Plot("He3gg-eta-theta-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr,"ET3"))
            .Hist(Hist(MC, r, histpath_central_reconstr,"ET4"))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'";
            Plot("He3gg-ggmm-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM4"), "before cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM5"), "after cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [2.2:3.4]"
                    << "set xlabel '2gamma missing mass, GeV'";
            Plot("He3gg-ggim-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "GIM5"), "before cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GIM6"), "after cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass + Q, GeV'"<< "set xrange [0:0.8]";
            Plot("He3gg-tim-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "TIM6-AllBins"))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            PlotHist2d(sp2, "He3gg-he3mm-he3me-mc" + r)
                .Distr(Hist2d(MC,r,histpath_central_reconstr,"He3MME6"))
                   << "set key on"
                   << "set xlabel '3He missing mass, GeV'"
                   << "set ylabel '3He  missing energy, GeV'";
    }
    {
            Plot("He3gg-he3mm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM0"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM1"), "cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM2"), "2 gammas required")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'"<< "set xrange [0.45:0.57]";

            PlotHist2d(sp2, "He3gg-alpha1-data")
                .Distr(Hist2d(DATA,"All",histpath_central_reconstr,"GGangle2").Scale(5,5))
                   << "set key on"
                   << "set xlabel 'sin alpha'"
                   << "set ylabel 'cos alpha'";
            PlotHist2d(sp2, "He3gg-alpha2-data")
                .Distr(Hist2d(DATA,"All",histpath_central_reconstr,"GGangle3").Scale(5,5))
                   << "set key on"
                   << "set xlabel 'sin alpha'"
                   << "set ylabel 'cos alpha'";

            Plot("He3gg-eta-theta-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"ET3"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"ET4"))
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'";
            Plot("He3gg-ggmm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM4"), "before cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM5"), "after cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'"<< "set xrange [2.2:3.4]";
            Plot("He3gg-ggim-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM5"), "before cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM6"), "after cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass + Q, GeV'"<< "set xrange [0:0.8]";

            const auto TIM=Hist(DATA, "All", histpath_central_reconstr, "TIM6-AllBins");
            Plot("He3gg-tim-data")
            .Hist(TIM, "IM and MM cuts")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                     << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";

            Plot("He3gg-dt-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt2"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He3gg-t-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t2"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";
            Plot("He3gg-dt-data-end")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt6"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He3gg-t-data-end")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t6"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";

            for(size_t a_t=0;a_t<suffix.size();a_t++){
                const double cutpos=-0.06+0.02*a_t;
                const double high=TIM.TransponateAndSort().right().X().max();
                Plot("He3gg-tim-data"+suffix[a_t]).Hist(TIM)
                    .Line({make_point(cutpos,0.),make_point(cutpos,high)})
                        << "set key on" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                        << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            }
            PlotHist2d(sp2, "He3gg-he3mm-he3me-data")
                .Distr(Hist2d(DATA,"All",histpath_central_reconstr,"He3MME6"))
                   << "set key on"
                   << "set xlabel '3He missing mass, GeV'"
                   << "set ylabel '3He  missing energy, GeV'";
    }

    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
            << "Q in [" << setprecision(3)<< Q.min() << "; " << Q.max() << "] MeV").str();
        cout<<Qmsg << " plots"<<endl;
        const auto TIM=Hist(DATA, "All", histpath_central_reconstr, string("TIM6-Bin-") + to_string(bin_num));
        for(size_t a_t=0;a_t<suffix.size();a_t++){
            const double cutpos=-0.06+0.02*a_t;
            const auto TIM_c=TIM.XRange(cutpos,2);
            cout<<Qmsg<< ";"<<suffix[a_t] <<endl;
            for(size_t i = 0; i < reaction.size(); i++){ 
                const auto &r = reaction[i];
                cout<<Qmsg << " acceptance "<<r<<endl;
                const auto MC_TIM=Hist(MC, r, histpath_central_reconstr, "TIM6-Bin-"+to_string(bin_num));
                hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
                const auto &N = Norm[bin_num].Y();
                if (N.Above(0)) acc[a_t][i] << make_point(Q, std_error(MC_TIM.XRange(cutpos,2).TotalSum().val())/N);
                else acc[a_t][i] << make_point(Q, 0.0);
            }
            cout<<Qmsg<< ";"<<suffix[a_t] << "; events count"<<endl;
            ev_am[a_t]<<make_point(Q,std_error(TIM_c.TotalSum().val()));
        }

    }
    const auto luminosity = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc"));
    const auto luminosity_he = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYf"));
    const auto true_he3eta = luminosity_he
        *extend_hist<2,2>(hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed")))
            .XRange(luminosity_he.left().X().min(),luminosity_he.right().X().max())
        /trigger_he3_forward.scaling;
    const double branching_ratio=0.39;

    for(size_t a_t=0;a_t<suffix.size();a_t++){
        cout<<suffix[a_t]<< " saving"<<endl;
        Plot accplot("He3gg-acceptance-final"+suffix[a_t]);
        accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:0.3]" << "set xrange [-70:30]"
            << "set key on";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto ac = acc[a_t][i].YRange(0.0000001, INFINITY);
            if (ac.size() > 0) {
                accplot.Hist(ac, reaction[i]);
                Plot("He3gg-acceptance-final"+suffix[a_t]+"-"+reaction[i])
                .Hist(acc[a_t][i],"","He3gg-acceptance"+suffix[a_t]+"-"+to_string(i))
                    << "set title 'Acceptance "+reaction[i]+"'"
                    << "set xlabel 'Q, MeV'"
                    << "set ylabel 'Acceptance, percents'"
                    << "set xrange [-70:30]";
            }
        }
        cout<<suffix[a_t]<< " fitting"<<endl;
        const auto known_events = (true_he3eta*branching_ratio)
            *extend_hist<2,2>(acc[a_t][3]).XRange(true_he3eta.left().X().min(),true_he3eta.right().X().max());
        Plot("He3gg-events-final"+suffix[a_t])
            .Hist(ev_am[a_t],"data","He3gg-data"+suffix[a_t])
            .Hist_2bars<1,2>(known_events,"3He+eta")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
        Plot("He3gg-events-final"+suffix[a_t]+"-bound")
            .Hist(ev_am[a_t].XRange(-70,2.5),"data, below threshold")
            .Hist_2bars<1,2>(extend_hist<1,2>(ev_am[a_t]).XRange(12.5,30)-known_events,"data-3Heeta, upper threshold")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]";
        const auto data_shape=(
            extend_hist<1,2>(ev_am[a_t])*trigger_he3_forward.scaling/(extend_hist<2,2>(acc[a_t][2])*luminosity)
        ).XRange(-70,2.5);
        Plot("He3gg-events-norm"+suffix[a_t]+"-bound-nofit")
            .Hist_2bars<1,2>(data_shape,"Data")
                << "set xlabel 'Q, MeV'" << "set key on"
                << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]";
    }
    cout<<"Final plots"<<endl;
    Plot("He3gg-tube-acc").Hist(
        Hist(MC, reaction[0], histpath_central_reconstr, "Events0")
        /Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference")
    )
            << "set xlabel 'Q, MeV'" << "set xrange [-70:30]"
            << "set ylabel 'Acceptance, n.d.'" << "set yrange [0:1]"
            << "set title 'How many helium ions from mesic nuclei decay would be detected'";
    cout<<"END"<<endl;
}
