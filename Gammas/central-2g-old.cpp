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
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma-old");
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound1-2g","bound2-2g","bound3-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("All");
    const hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";

    hist<> ev_am;
    vector<hist<>> acc;
    for (size_t j = 0; j < reaction.size(); j++)
        acc.push_back(hist<>());

    for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            cout<<"All-bins MC plots "<<r<<endl;
            Plot("He3gg-old-he3mm-mc-" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr,"old_He3MM0"))
            .Hist(Hist(MC, r, histpath_central_reconstr,"old_He3MM1"), "cut")
            .Hist(Hist(MC, r, histpath_central_reconstr,"old_He3MM2"), "2 gammas required")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [0.45:0.57]"
                    << "set xlabel '3He missing mass - Q, GeV'";
            Plot("He3gg-old-ggmm-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "old_GMM2"), "before cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, "old_GMM3"), "after cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [2.2:3.4]"
                    << "set xlabel '2gamma missing mass, GeV'";
            Plot("He3gg-old-ggim-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "old_GIM3"))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'"<< "set xrange [0:0.8]";
            Plot("He3gg-old-tim-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "old_TIM3"))
            .Hist(Hist(MC, r, histpath_central_reconstr, "old_TIM4"))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
    }
    {
            Plot("He3gg-old-he3mm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_He3MM0"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_He3MM1"), "cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_He3MM2"), "2 gammas required")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'"<< "set xrange [0.45:0.57]";
            Plot("He3gg-old-ggmm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_GMM2"), "before cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_GMM3"), "after cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'"<< "set xrange [2.2:3.4]";
            Plot("He3gg-old-ggim-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_GIM3"))
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'"<< "set xrange [0:0.8]";
            Plot("He3gg-old-tim-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_TIM3"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_TIM4"))
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                     << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            Plot("He3gg-old-dt-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_dt4"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_dt5"),"cut")
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He3gg-old-t-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "old_t5-AllBins"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";
    }

    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
            << "Q in [" << setprecision(3)<< Q.min() << "; " << Q.max() << "] MeV").str();
        cout<<Qmsg << " plots"<<endl;
        const auto DataT=Hist(DATA, "All",histpath_central_reconstr, string("old_t5-Bin-") + to_string(bin_num)).Scale(10);
        const auto DataTCut=DataT.XRange(-10,30);
            for(size_t i = 0; i < reaction.size(); i++){ 
                const auto &r = reaction[i];
                cout<<Qmsg << " acceptance "<<r<<endl;
                const auto MCT=Hist(MC, r, histpath_central_reconstr, "old_t5-Bin-"+to_string(bin_num));
                hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
                const auto &N = Norm[bin_num].Y();
                if (N.Above(0)) acc[i] << make_point(Q, std_error(MCT.TotalSum().val())/N);
                else acc[i] << make_point(Q, 0.0);
            }
        Plot(Q.Contains(-11) ? "He3gg-old-t-data-part" : "")
            .Hist(DataT).Hist(DataTCut)
            << "set title 'Data " + runmsg +" " + Qmsg + "'"  << "set yrange [0:]"
            << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";

        ev_am<<make_point(Q,std_error(DataTCut.TotalSum().val()));
    }


    Plot accplot("He3gg-old-acceptance-final-20");
    accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:0.2]" << "set xrange [-70:30]"
            << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
            const auto ac = acc[i].YRange(0.0000001, INFINITY);
            if (ac.size() > 0) {
                accplot.Hist(ac, reaction[i]);
                Plot("He3gg-old-acceptance-final-20-"+reaction[i])
                .Hist(acc[i],"","He3gg-acceptance-20-"+to_string(i))
                    << "set title 'Acceptance "+reaction[i]+"'"
                    << "set xlabel 'Q, MeV'"
                    << "set ylabel 'Acceptance, percents'"
                    << "set xrange [-70:30]";
            }
    }
    const ext_hist<2> luminosity_he = Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYf");
    const auto true_he3eta = luminosity_he
        *extend_hist<2,2>(hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed")))
            .XRange(luminosity_he.left().X().min(),luminosity_he.right().X().max())
        /trigger_he3_forward.scaling;
    const double branching_ratio=0.39;
    const auto known_events = (true_he3eta*branching_ratio)
        *extend_hist<2,2>(acc[3]).XRange(true_he3eta.left().X().min(),true_he3eta.right().X().max());

    Plot("He3gg-old-events-final-20")
            .Hist(ev_am,"data").Hist_2bars<1,2>(known_events,"3He+eta")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
    const ext_hist<2> luminosity = Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc");
    const auto ev_norm=(
        extend_hist<1,2>(ev_am)*trigger_he3_forward.scaling/(extend_hist<2,2>(acc[2])*luminosity)
    ).XRange(-40,2.5);

    Plot("He3gg-old-events-final-20-part")
            .Hist_2bars<1,2>(ev_norm)
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-40:2.5]"
                << "set ylabel 'events normalized, nb'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
    cout<<"END"<<endl;
}
