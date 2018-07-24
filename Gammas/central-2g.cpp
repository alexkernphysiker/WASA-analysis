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
#include <Parameters/parameters.h>
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma");
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound1-2g","bound2-2g","bound3-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
    const hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const auto runs = PresentRuns("All");
    const string runmsg = (runs.first==runs.second)?"":"("+to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs)";
    cout<<"First plots"<<endl;
    Plot gE("He3gg-gamma-energy-theory",3), gEc("He3gg-gamma-energy-theory-cut",3),
         gEd("He3gg-gamma-energy-data",3);
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
    Plot("He3gg-gamma-count",3).Hist(Hist(DATA, "All", histpath_central_reconstr, "GammaCount"));

    hist<> ev_am;
    vector<hist<>> acceptance;
    vector<hist<>> acc;

    for (size_t j = 0; j < reaction.size(); j++)
        acc.push_back(hist<>());

    for (size_t i = 0; i < reaction.size(); i++) {
            const int N=10000000;
            const auto &r = reaction[i];
            cout<<"All-bins MC plots "<<r<<endl;
            const auto he3mm0=Hist(MC, r, histpath_central_reconstr,"He3MM0")/N;
            Plot("He3gg-he3mm-mc-" + r,5)
            .Hist(he3mm0).Hist(Hist(MC, r, histpath_central_reconstr,"He3MM1")/N)
            .Line(Points<>{{getParameter(he3mm_cut),0.0},{getParameter(he3mm_cut),he3mm0.TransponateAndSort().right().X().max()*1.5}},"cut")
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM2")/N, "2 gammas required")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [0.45:0.57]"
                    << "set xlabel '3He missing mass - Q, GeV'"<<"set ylabel 'Efficiensy, a.u.'";

            const auto ggcos=Hist(MC, r, histpath_central_reconstr,"GGcos0")/N;
            Plot("He3gg-cos-mc" + r,5)
            .Hist(ggcos).Hist(Hist(MC, r, histpath_central_reconstr,"GGcos3")/N)
            .Line(Points<>{{0.0,0.0},{0.0,hist<>(ggcos.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(alpha)'"<<"set ylabel 'Efficiensy, a.u.'";
            const auto eta_theta=Hist(MC, r, histpath_central_reconstr,"ET3")/N;
            Plot("He3gg-eta-theta-mc" + r,5)
            .Hist(eta_theta).Hist(Hist(MC, r, histpath_central_reconstr,"ET4")/N)
            .Line(Points<>{{getParameter(eta_theta_thr),0.0},{getParameter(eta_theta_thr),hist<>(eta_theta.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'"<<"set ylabel 'Efficiensy, a.u.'";
            const auto ggmm=Hist(MC, r, histpath_central_reconstr, "GMM4")/N;
            Plot("He3gg-ggmm-mc" + r,5)
            .Hist(ggmm).Hist(Hist(MC, r, histpath_central_reconstr, "GMM5")/N)
            .Line(Points<>{{getParameter(gamma_mm_lo),hist<>(ggmm.Transponate()).right().X().max()*1.5},{getParameter(gamma_mm_lo),0.0},
                           {getParameter(gamma_mm_hi),0.0},{getParameter(gamma_mm_hi),hist<>(ggmm.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [2.2:3.4]"
                    << "set xlabel '2gamma missing mass, GeV'"<<"set ylabel 'Efficiensy, a.u.'";
            const auto ggim=Hist(MC, r, histpath_central_reconstr, "GIM5")/N;
            Plot("He3gg-ggim-mc" + r,5)
            .Hist(ggim).Hist(Hist(MC, r, histpath_central_reconstr, "GIM6")/N)
            .Line(Points<>{{getParameter(gamma_im_lo),hist<>(ggim.Transponate()).right().X().max()*1.5},{getParameter(gamma_im_lo),0.0},
                           {getParameter(gamma_im_hi),0.0},{getParameter(gamma_im_hi),hist<>(ggim.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<<"set ylabel 'Efficiensy, a.u.'"
                    << "set xlabel '2gamma invariant mass - Q, GeV'"<< "set xrange [0:0.8]";
            Plot("He3gg-tim-mc" + r,5)
            .Hist(Hist(MC, r, histpath_central_reconstr, "TIM6-AllBins")/N)
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'"<<"set ylabel 'Efficiensy, a.u.'";
    }
    {
            const auto he3mm0=Hist(DATA, "All", histpath_central_reconstr,"He3MM0");
            Plot("He3gg-he3mm-data",5)
            .Hist(he3mm0).Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM1"))
            .Line(Points<>{{getParameter(he3mm_cut),0.0},{getParameter(he3mm_cut),hist<>(he3mm0.Transponate()).right().X().max()*1.5}},"cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM2"), "2 gammas required")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel '3He missing mass - Q, GeV'"<< "set xrange [0.45:0.57]";

            const auto ggcos=Hist(DATA, "All", histpath_central_reconstr,"GGcos0");
            Plot("He3gg-cos-data",5)
            .Hist(ggcos).Hist(Hist(DATA, "All", histpath_central_reconstr,"GGcos3"))
            .Line(Points<>{{0.0,0.0},{0.0,hist<>(ggcos.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title 'Data"+runmsg+"'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(alpha)'"<<"set ylabel 'Events, n.d.'";
            const auto eta_theta=Hist(DATA,"All", histpath_central_reconstr,"ET3");
            Plot("He3gg-eta-theta-data",5)
            .Hist(eta_theta).Hist(Hist(DATA, "All", histpath_central_reconstr,"ET4"))
            .Line(Points<>{{getParameter(eta_theta_thr),0.0},{getParameter(eta_theta_thr),hist<>(eta_theta.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'"<<"set ylabel 'Events, n.d.'";
            const auto ggmm=Hist(DATA,"All", histpath_central_reconstr, "GMM4");
            Plot("He3gg-ggmm-data",5)
            .Hist(ggmm).Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM5"))
            .Line(Points<>{{getParameter(gamma_mm_lo),hist<>(ggmm.Transponate()).right().X().max()*1.5},{getParameter(gamma_mm_lo),0.0},
                           {getParameter(gamma_mm_hi),0.0},{getParameter(gamma_mm_hi),hist<>(ggmm.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'"<< "set xrange [2.2:3.4]"<<"set ylabel 'Events, n.d.'";
            const auto ggim=Hist(DATA,"All", histpath_central_reconstr, "GIM5");
            Plot("He3gg-ggim-data",5)
            .Hist(ggim).Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM6"))
            .Line(Points<>{{getParameter(gamma_im_lo),hist<>(ggim.Transponate()).right().X().max()*1.5},{getParameter(gamma_im_lo),0.0},
                           {getParameter(gamma_im_hi),0.0},{getParameter(gamma_im_hi),hist<>(ggim.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel '2gamma invariant mass - Q, GeV'"<< "set xrange [0:0.8]";

            Plot("He3gg-tim-data",5)
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "TIM6-AllBins"))
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-0.3:0.3]"
                     << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'"<<"set ylabel 'Events, n.d.'";

            const auto DT=Hist(DATA, "All", histpath_central_reconstr, "dt00");
            Plot("He3gg-dt-data",5)
            .Hist(DT).Hist(Hist(DATA, "All", histpath_central_reconstr, "dt0"))
            .Line(Points<>{{getParameter(time_dt),0.0},{getParameter(time_dt),hist<>(DT.Transponate()).right().X().max()*1.5}},"cut")
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            const auto T=Hist(DATA, "All", histpath_central_reconstr, "t00");
            Plot("He3gg-t-data",5)
            .Hist(T).Hist(Hist(DATA, "All", histpath_central_reconstr, "t0"))
            .Line(Points<>{{getParameter(time_t1),hist<>(T.Transponate()).right().X().max()*1.5},{getParameter(time_t1),0.0},
                           {getParameter(time_t2),0.0},{getParameter(time_t2),hist<>(T.Transponate()).right().X().max()*1.5}},"cut")
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";

            Plot("He3gg-dt-data-end",5)
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt6"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on"<<"set ylabel 'Events, n.d.'";
            Plot("He3gg-t-data-end",5)
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t6"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on"<<"set ylabel 'Events, n.d.'";
            Plot("He3gg-cos-data-end",5)
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"GGcos6"),"data")
                    << "set key on" << "set title 'Data"+runmsg+"'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(alpha)'"<<"set ylabel 'Events, n.d.'";
            Plot("He3gg-eta-theta-data-end",5)
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"ET6"),"data")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'"<<"set ylabel 'Events, n.d.'";
    }

    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
            << "Q in [" << setprecision(3)<< Q.min() << "; " << Q.max() << "] MeV").str();
        cout<<Qmsg << " plots"<<endl;
        const auto TIM=Hist(DATA, "All", histpath_central_reconstr, string("TIM6-Bin-") + to_string(bin_num)).XRange(-0.3,0.3);
        for(size_t i = 0; i < reaction.size(); i++){ 
            const auto &r = reaction[i];
            cout<<Qmsg << " acceptance "<<r<<endl;
            const auto MC_TIM=Hist(MC, r, histpath_central_reconstr, "TIM6-Bin-"+to_string(bin_num)).XRange(-0.3,0.3);
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

        Plot accplot("He3gg-acceptance-final",5),accplot2("He3gg-acceptance-bg",5);
        accplot << "set title 'Efficiency'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Efficiency, n.d.'"
            << "set yrange [0:0.4]" << "set xrange [-70:30]"
            << "set key on";
        accplot2 << "set title 'Efficiency'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Efficiency, n.d.'"
            << "set yrange [0:0.005]" << "set xrange [-70:30]"
            << "set key on";
        for (size_t i = 0; i < reaction.size(); i++) {
            accplot.Hist(acc[i], reaction[i]);
        }
        for (size_t i = 4; i < reaction.size(); i++) {
            accplot2.Hist(acc[i], reaction[i]);
        }
        const auto known_events = (true_he3eta*branching_ratio)*extend_hist<2,2>(acc[3]);
        Plot("He3gg-events-final",5)
            .Hist(ev_am,"data")
            .Hist_2bars<1,2>(known_events,"3He+eta")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'Events, n.d.'" << "set yrange [0:]"
                << "set title '"+runmsg+"'"<<"set key left top";
        Plot("He3gg-events-final-cut",5)
            .Hist(ev_am,"data, below threshold")
            .Hist_2bars<1,2>(extend_hist<1,2>(ev_am)-known_events,"data-3Heeta")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'Events, n.d.'" << "set yrange [0:]";
        const auto data_shape=(extend_hist<1,2>(ev_am)-known_events)*trigger_he3_forward.scaling/(extend_hist<2,2>(acc[0])*luminosity); 
        Plot("He3gg-events-norm",5)
            .Hist_2bars<1,2>(data_shape.XRange(-70,10),"Data statistical", "Data systematic")
                << "set xlabel 'Q, MeV'" << "set key on"<< "set title '"+runmsg+"'"<< "set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:]";
        Plot("He3gg-events-norm-light",5)
            .Hist(wrap_hist(data_shape).XRange(-70,10))
                << "set xlabel 'Q, MeV'" << "set key on"<< "set title '3He+2gamma"+runmsg+"'"<< "set xrange [-70:10]"
                << "set ylabel 'Normalized events, nb'" << "set yrange [0:]";

    cout<<"Final plots"<<endl;
    Plot("He3gg-tube-acc",5).Hist(
        Hist(MC, reaction[0], histpath_central_reconstr, "Events0")
        /Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference")
    )
            << "set xlabel 'Q, MeV'" << "set xrange [-70:30]"
            << "set ylabel 'Efficiency, n.d.'" << "set yrange [0:1]"
            << "set title 'How many helium ions from mesic nuclei decay would be detected'";
    cout<<"END"<<endl;
}
