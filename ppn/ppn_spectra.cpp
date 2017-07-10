// this file is distributed under 
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
    const auto runs=PresentRuns("L");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";
    const string th1="'Theta_1, deg'",th2="'Theta_2, deg'",e1="'Edep_1, GeV'",e2="'Edep_2, GeV'",
    thth="'Theta_1+1.6Theta_2, deg'",planarity="'|Phi_1-Phi_2-180^o|, deg'";
    const hist<> norm=Hist(MC,"ppn_qf",{"Histograms","elastic"},"0-Reference");
    const hist<> norm_pd=Hist(MC,"pd",{"Histograms","elastic"},"0-Reference");
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"ppn");
    Plot<>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_0")/norm_pd.TotalSum().val(),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_0")/norm.TotalSum().val(),"ppn_{sp}")<<"set key on"
    <<"set title 'Planarity. MC'"<<"set yrange [0:]"<<"set xlabel "+planarity;
    Plot<>()
    .Hist(Hist(MC,"pd",{"Histograms","elastic"},"pair_phi_diff_1")/norm_pd.TotalSum().val(),"pd")
    .Hist(Hist(MC,"ppn_qf",{"Histograms","elastic"},"pair_phi_diff_1")/norm.TotalSum().val(),"ppn_{sp}")<<"set key on"
    <<"set title 'Planarity. MC. Cut'"<<"set yrange [0:]"<<"set xlabel "+planarity;
    Plot<>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"pair_phi_diff_0"))
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"pair_phi_diff_1"))
    <<"set title 'Planarity. Data "+runmsg+"'"<<"set yrange [0:]"<<"set xlabel "+planarity;

    PlotHist2d<>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'MC pd'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'MC ppn_{sp}'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_1"))
    <<"set zrange [0:]"<<"set title 'Data "+runmsg+"'"<<"set xlabel "+th1<<"set ylabel "+th2;


    PlotHist2d<>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC pd'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC ppn_{sp}'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_t_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. Data "+runmsg+"'"<<"set xlabel "+th1<<"set ylabel "+th2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"t_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC pd'"<<"set xlabel "+th1<<"set ylabel "+e1;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"t_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1 MC ppn_{sp}'"<<"set xlabel "+th1<<"set ylabel "+e1;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"t_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. Data "+runmsg+"'"<<"set xlabel "+th1<<"set ylabel "+e1;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"pd",{"Histograms","elastic"},"e_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC pd'"<<"set xlabel "+e1<<"set ylabel "+e2;
    PlotHist2d<>(sp2).Distr(Hist2d(MC,"ppn_qf",{"Histograms","elastic"},"e_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. MC ppn_{sp}'"<<"set xlabel "+e1<<"set ylabel "+e2;
    PlotHist2d<>(sp2).Distr(Hist2d(DATA,"L",{"Histograms","elastic"},"e_vs_e_22"))
    <<"set zrange [0:]"<<"set title 'Cut1. Data "+runmsg+"'"<<"set xlabel "+e1<<"set ylabel "+e2;
    Plot<>()
    .Line(Hist(MC,"pd",{"Histograms","elastic"},"theta_sum_22-AllBins").toLine()/norm_pd.TotalSum().val(),"pd")
    .Line(Hist(MC,"ppn_qf",{"Histograms","elastic"},"theta_sum_22-AllBins").toLine()/norm.TotalSum().val(),"ppn_{sp}")
    <<"set title 'MC'"<<"set key on"<<"set yrange [0:]"<<"set xrange [0:300]"<<"set xlabel "+thth;
    Plot<>()
    .Hist(Hist(DATA,"L",{"Histograms","elastic"},"theta_sum_22-AllBins"))
    <<"set title 'Data "+runmsg+"'"<<"set key on"<<"set yrange [0:]"<<"set xrange [0:300]"<<"set xlabel "+thth;



    hist<> acceptance,acceptance_pd;
    for(size_t bin_num=0,bin_count=norm.size();bin_num<bin_count;bin_num++){
	const auto&Q=norm[bin_num].X();
	const auto&N=norm[bin_num].Y();
	const auto&N_pd=norm_pd[bin_num].Y();
	const string Qmsg=static_cast<stringstream&>(stringstream()
	    <<"Q in ["<<setprecision(3)
	    <<Q.min()<<"; "<<Q.max()<<"] MeV"
	).str();

	const hist<> mc_ppn=Hist(MC,"ppn_qf",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num));
	const hist<> mc_pd=Hist(MC,"pd",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num));
	const hist<> nmc_ppn=mc_ppn/(N*mc_ppn[0].X().uncertainty()*2.);
	const hist<> nmc_pd=mc_pd/(N_pd*mc_pd[0].X().uncertainty()*2.);
	acceptance<<point<value<>>(Q,mc_ppn.TotalSum()/N);
	acceptance_pd<<point<value<>>(Q,mc_pd.TotalSum()/N_pd);
	const hist<> data=Hist(DATA,"L",{"Histograms","elastic"},string("theta_sum_22-Bin-")+to_string(bin_num));

	Plot<>().Hist(nmc_ppn,"ppn_{sp}").Hist(nmc_pd,"pd")
	<<"set key on"<<"set title 'MC "+Qmsg+"'"<<"set yrange [0:0.008]"
	<<"set xlabel "+thth<<"set ylabel 'Differential acceptance, deg^{-1}'";
	Plot<>()
	.Hist(data,"data")
	<<"set key on"<<"set title 'Data "+Qmsg+" "+runmsg+"'"<<"set yrange [0:]"
	<<"set xlabel "+thth<<"set ylabel 'Count'";


    }
    Plot<>().Hist(acceptance,"ppn_{sp}").Hist(acceptance_pd,"pd")<<"set key on"
    <<"set title 'Acceptance'"<<"set yrange [0:]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
}
