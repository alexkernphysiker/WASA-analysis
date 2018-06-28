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
    vector<string> reaction = {"bound1-2g","bound1-2g","bound1-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
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

    hist<> ev_am;
    vector<hist<>> acceptance;
    vector<hist<>> acc;
    const double cut_pos=-0.200;

    for (size_t j = 0; j < reaction.size(); j++)
        acc.push_back(hist<>());

    for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            cout<<"All-bins MC plots "<<r<<endl;
            Plot("He3gg-he3mm-mc-" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM0"))
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM1"), "cut")
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM2"), "2 gammas required")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [0.45:0.57]"
                    << "set xlabel '3He missing mass - Q, GeV'";

            Plot("He3gg-cos-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr,"GGcos0"))
            .Hist(Hist(MC, r, histpath_central_reconstr,"GGcos3"))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(alpha)'";
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
            .Hist(Hist(MC, r, histpath_central_reconstr, "TIM6-AllBins").XRange(cut_pos,0.2))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            Plot("He3gg-he3me-mc"+r)
                .Hist(Hist(MC, r, histpath_central_reconstr, "He3ME6-AllBins"))
                   << "set key on"<<"set yrange [0:]"<< "set title '"+r+"'"
                   << "set xlabel '3He missing energy, GeV'";
    }
    {
            Plot("He3gg-he3mm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM0"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM1"), "cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM2"), "2 gammas required")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'"<< "set xrange [0.45:0.57]";

            Plot("He3gg-cos-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"GGcos0"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"GGcos3"))
                    << "set key on" << "set title 'Data"+runmsg+"'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(alpha)'";
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

            Plot("He3gg-tim-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "TIM6-AllBins"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "TIM6-AllBins").XRange(cut_pos,0.2))
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                     << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";

            Plot("He3gg-dt-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt00"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt0"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He3gg-t-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t00"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t0"))
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
            Plot("He3gg-cos-data-end")
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"GGcos6"),"data")
                    << "set key on" << "set title 'Data"+runmsg+"'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(alpha)'";
            Plot("He3gg-eta-theta-data-end")
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"ET6"),"data")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'";

            Plot("He3gg-he3me-data")
                .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3ME6-AllBins"))
                   << "set key on"<<"set yrange [0:]"<< "set title 'Data " + runmsg + "'"
                   << "set xlabel '3He missing energy, GeV'";
    }

    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
            << "Q in [" << setprecision(3)<< Q.min() << "; " << Q.max() << "] MeV").str();
        cout<<Qmsg << " plots"<<endl;
        const auto TIM=Hist(DATA, "All", histpath_central_reconstr, string("TIM6-Bin-") + to_string(bin_num)).XRange(cut_pos,0.2);
        for(size_t i = 0; i < reaction.size(); i++){ 
            const auto &r = reaction[i];
            cout<<Qmsg << " acceptance "<<r<<endl;
            const auto MC_TIM=Hist(MC, r, histpath_central_reconstr, "TIM6-Bin-"+to_string(bin_num)).XRange(cut_pos,0.2);
            hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) acc[i] << make_point(Q, std_error(MC_TIM.TotalSum().val())/N);
            else acc[i] << make_point(Q, 0.0);
        }
        ev_am<<make_point(Q,std_error(TIM.TotalSum().val()));

    }
    const auto luminosity = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc"));
    const auto luminosity_he = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc")).XRange(10,30);
    auto true_he3eta = luminosity_he
        *extend_hist<2,2>(hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed")).XRange(10,30))/trigger_he3_forward.scaling;
    while(true_he3eta.left().X().min()>-70.)true_he3eta<<make_point(value<>(true_he3eta.left().X().val()-2.5,2.5),0);
    const auto branching_ratio=uncertainties(0.393,0,0.003);

        Plot accplot("He3gg-acceptance-final"),accplot2("He3gg-acceptance-bg");
        accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:0.4]" << "set xrange [-70:30]"
            << "set key on";
        accplot2 << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:0.005]" << "set xrange [-70:30]"
            << "set key on";
        for (size_t i = 0; i < reaction.size(); i++) {
            accplot.Hist(acc[i], reaction[i]);
        }
        for (size_t i = 4; i < reaction.size(); i++) {
            accplot2.Hist(acc[i], reaction[i]);
        }
        const auto known_events = (true_he3eta*branching_ratio)*extend_hist<2,2>(acc[3]);
        Plot("He3gg-events-final")
            .Hist(ev_am,"data")
            .Hist_2bars<1,2>(known_events,"3He+eta")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]"
                << "set title '"+runmsg+"'"<<"set key left top";
        Plot("He3gg-events-final-cut")
            .Hist(ev_am,"data, below threshold")
            .Hist_2bars<1,2>(extend_hist<1,2>(ev_am)-known_events,"data-3Heeta")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]";
        const auto data_shape=(extend_hist<1,2>(ev_am)-known_events)*trigger_he3_forward.scaling/(extend_hist<2,2>(acc[0])*luminosity); 
        Plot("He3gg-events-norm")
            .Hist_2bars<1,2>(data_shape.XRange(-70,10),"Data statistical", "Data systematic")
                << "set xlabel 'Q, MeV'" << "set key on"<< "set title '"+runmsg+"'"<< "set xrange [-70:10]"
                << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]";

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
