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
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-6gamma-with-3he");
    const vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas6"};
    const vector<string> reaction = {"bound1-6g","bound2-6g","bound3-6g", "He3eta-6g", "He3pi0pi0pi0"};
    const vector<string> rname = {"Bound state","Bound state","Bound state", "pd->^3Heη", "pd->^3He3π^0"};
    hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const auto runs = PresentRuns("All");
    const string runmsg = (runs.first==runs.second)?"":"("+to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs)";

    const int N=10000000;
    
    Plot gE("He36g-gamma-energy-theory",4), gEc("He36g-gamma-energy-theory-cut",4),
         gEd("He36g-gamma-energy-data",4);
    gE << "set key on" << "set title 'γ energy. MC'"
       << "set yrange [0:]" << "set xlabel 'γ energy, GeV'";
    gEc << "set key on" << "set title 'γ energy. MC. Cut'"
        << "set yrange [0:]" << "set xlabel 'γ energy, GeV'";
    gEd << "set key on" << "set title 'γ energy. Data'"
        << "set yrange [0:]" << "set xlabel 'γ energy, GeV'";
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

    {
        Plot theory("He36g-IMPiDiff-mc",4),experiment("He36g-IMPiDiff-data",4);
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            const auto &rn = rname[i];
	    const auto imdiff=Hist(MC, r, histpath_central_reconstr, "GMMPDiff4")/N;
            if (i == 1) {
                Plot("He36g-IMPiDiff-bound-mc",3)
                .Hist(imdiff).Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff5")/N)
                .Line(Points<>{{getParameter(three_pi0),0.0},{getParameter(three_pi0),hist<>(imdiff.Transponate()).right().X().max()*1.5}})
                    << "set key on" << "set yrange [0:]"<<"set xlabel 'D, GeV^2/c^4'";
            }
            if((i!=0)&&(i!=2))theory.Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff5")/N, rn)
            .Line(Points<>{{getParameter(three_pi0),0.0},{getParameter(three_pi0),hist<>(imdiff.Transponate()).right().X().max()*1.5}});
        }
        theory << "set key on" << "set yrange [0:]"<<"set xrange [0:0.15]"
	    <<"set ylabel 'Efficiency, a.u.'"<<"set title 'Monte Carlo'"<<"set xlabel 'D, GeV^2/c^4'";
        const auto imdiff=Hist(DATA, "All", histpath_central_reconstr, "GMMPDiff4");
        experiment
            .Hist(imdiff).Hist(Hist(DATA, "All", histpath_central_reconstr, "GMMPDiff5"))
            .Line(Points<>{{getParameter(three_pi0),0.0},{getParameter(three_pi0),hist<>(imdiff.Transponate()).right().X().max()*1.5}})
                << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
		<<"set xrange [0:0.15]"<<"set ylabel 'Events, n.d.'"<<"set xlabel 'D, GeV^2/c^4'";
    }
    for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            const auto &rn = rname[i];
            cout<<"All-bins MC plots "<<r<<endl;
            const auto cosb=Hist(MC, r, histpath_central_reconstr,"cosi2")
                        +Hist(MC, r, histpath_central_reconstr,"cosj2")
                        +Hist(MC, r, histpath_central_reconstr,"cosk2");
            const auto cosa=Hist(MC, r, histpath_central_reconstr,"cosi3")
                        +Hist(MC, r, histpath_central_reconstr,"cosj3")
                        +Hist(MC, r, histpath_central_reconstr,"cosk3");
            Plot("He36g-cos-mc"+r,5)
            .Hist(cosb).Hist(cosa)
                    << "set key on" << "set title '" + r + "'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(γ-γ)'";
            const auto he3mm0=Hist(MC, r, histpath_central_reconstr,"He3MM0")/N;
            Plot("He36g-he3mm-mc-" + r,4)
            .Hist(he3mm0).Hist(Hist(MC, r, histpath_central_reconstr,"He3MM1")/N)
            .Line(Points<>{{getParameter(he3mm_cut),0.0},{getParameter(he3mm_cut),hist<>(he3mm0.Transponate()).right().X().max()*1.5}},"cut")
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM2")/N, "6γ required")
                    << "set key on" << "set title '"+rn+"'" << "set yrange [0:]"<<"set ylabel 'Efficiensy, a.u.'"
                    << "set xlabel '3He missing mass - Q_{3Heη}, GeV/c^2'"<< "set xrange [0.45:0.57]";
            const auto eta_theta=Hist(MC, r, histpath_central_reconstr,"ET3")/N;
            Plot("He36g-eta-theta-mc"+r,4)
            .Hist(eta_theta).Hist(Hist(MC, r, histpath_central_reconstr,"ET4")/N)
            .Line(Points<>{{getParameter(eta_theta_thr),0.0},{getParameter(eta_theta_thr),hist<>(eta_theta.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title '" + rn + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'θ(η) reconstructed, deg'"<<"set ylabel 'Efficiensy, a.u.'";
            const auto ggmm=Hist(MC, r, histpath_central_reconstr, "GMM5")/N;
            Plot("He36g-6gmm-mc" + r,4)
            .Hist(ggmm).Hist(Hist(MC, r, histpath_central_reconstr, "GMM6")/N)
            .Line(Points<>{{getParameter(gamma_mm_lo),hist<>(ggmm.Transponate()).right().X().max()*1.5},{getParameter(gamma_mm_lo),0.0},
                           {getParameter(gamma_mm_hi),0.0},{getParameter(gamma_mm_hi),hist<>(ggmm.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title '"+rn+"'" << "set yrange [0:]"<<"set ylabel 'Efficiensy, a.u.'"
                    << "set xlabel '6γ missing mass, GeV/c^2'"<< "set xrange [2.2:3.4]";
            const auto ggim=Hist(MC, r, histpath_central_reconstr, "GIM6")/N;
            Plot("He36g-6gim-mc" + r,4)
            .Hist(ggim).Hist(Hist(MC, r, histpath_central_reconstr, "GIM7")/N)
            .Line(Points<>{{getParameter(gamma_im_lo6),hist<>(ggim.Transponate()).right().X().max()*1.5},{getParameter(gamma_im_lo6),0.0},
                {getParameter(gamma_im_hi6),0.0},{getParameter(gamma_im_hi6),hist<>(ggim.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title '"+rn+"'" << "set yrange [0:]"<<"set ylabel 'Efficiensy, a.u.'"
                    << "set xlabel '6γ invariant mass - Q_{3Heη}, GeV/c^2'"<< "set xrange [0.0:1.0]";
            const auto lasthist=Hist(MC, r, histpath_central_reconstr, "TIM7-AllBins").Scale(4)/N;
            Plot("He36g-tim-mc" + r,4)
            .Hist(lasthist)
                    << "set key on" << "set title '"+rn+"'" << "set yrange [0:]"<< "set xrange [-0.2:0.15]"
                    << "set xlabel 'IM(^3He+6γ)-IM(p+d), GeV/c^2'"<<"set ylabel 'Efficiensy, a.u.'";
    }
    static Cache<string,hist<>> NORM;
    cout<<"========estimating background======="<<endl;
    const auto get_lh=[&histpath_central_reconstr](){
        auto lasthist=Hist(DATA, "All", histpath_central_reconstr, "TIM7-Bin-32");
        for(size_t bin_cnt=33;bin_cnt<40;bin_cnt++)
            lasthist+=Hist(DATA, "All", histpath_central_reconstr, "TIM7-Bin-"+to_string(bin_cnt));
        return lasthist.Scale(4);
    };
    const auto luminosity = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc"));
    const auto luminosity_he = luminosity.XRange(10,30);
    const auto branching_ratio=uncertainties(0.3257,0,0.0023);
    const auto he3eta_cs=hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed")).XRange(10,30);
    const auto get_lh_mc=[&histpath_central_reconstr,&reaction,&luminosity,&he3eta_cs,&branching_ratio](){
        const auto &N = NORM(reaction[3]+"_",[&histpath_central_reconstr,&reaction](){
          return Hist(MC, reaction[3], histpath_central_reconstr, "0-Reference","_");
        });
        hist<> lasthist;
        for(size_t bin_cnt=32;bin_cnt<40;bin_cnt++){
	    auto LH=
		(Hist(MC, reaction[3], histpath_central_reconstr, "TIM7-Bin-"+to_string(bin_cnt))/N[bin_cnt].Y())
		*(luminosity[bin_cnt].Y()*branching_ratio).wrap()
		*he3eta_cs[bin_cnt-32].Y()/trigger_he3_forward.scaling;
	    if(lasthist.size()==0)lasthist=LH;
	    else lasthist+=LH;
	}
        return lasthist.Scale(4);
    };
    const auto lasthist=get_lh();
    const auto he3eta_hist=get_lh_mc();
    const value<> S=lasthist.TotalSum()-he3eta_hist.TotalSum();
    Plot("He36g-tim-BG",5)
        .Hist(lasthist,"data")
	.Hist(he3eta_hist,"pd->^3Heη")
            << "set key on" << "set title 'Q є [10;30] MeV " + runmsg + "'" << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
            << "set xrange [-0.2:0.15]" << "set xlabel 'IM(^3He+6γ)-IM(p+d), GeV/c^2'";
    {
            const auto cosb=Hist(DATA, "All", histpath_central_reconstr,"cosi2")
                        +Hist(DATA, "All", histpath_central_reconstr,"cosj2")
                        +Hist(DATA, "All", histpath_central_reconstr,"cosk2");
            const auto cosa=Hist(DATA, "All", histpath_central_reconstr,"cosi3")
                        +Hist(DATA, "All", histpath_central_reconstr,"cosj3")
                        +Hist(DATA, "All", histpath_central_reconstr,"cosk3");
            Plot("He36g-cos-data",5)
            .Hist(cosb).Hist(cosa)
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(γ-γ), n.d.'";
            const auto he3mm0=Hist(DATA, "All", histpath_central_reconstr, "He3MM0")/1000.;
            Plot("He36g-he3mm-data",5)
            .Hist(he3mm0).Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM1")/1000.)
            .Line(Points<>{{getParameter(he3mm_cut),0.0},{getParameter(he3mm_cut),hist<>(he3mm0.Transponate()).right().X().max()*1.5}},"cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM2")/1000., "6γ required")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<<"set ylabel 'Events, 10^3'"
                    << "set xlabel '3He missing mass - Q_{3Heη}, GeV/c^2'"<< "set xrange [0.45:0.57]";
            const auto eta_theta=Hist(DATA, "All", histpath_central_reconstr,"ET3");
            Plot("He36g-eta-theta-data",5)
            .Hist(eta_theta).Hist(Hist(DATA, "All", histpath_central_reconstr,"ET4"))
            .Line(Points<>{{getParameter(eta_theta_thr),0.0},{getParameter(eta_theta_thr),hist<>(eta_theta.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'θ(η) reconstructed, deg'"<<"set ylabel 'Events, n.d.'";
            const auto ggmm=Hist(DATA, "All", histpath_central_reconstr, "GMM5");
            Plot("He36g-6gmm-data",5)
            .Hist(ggmm).Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM6"))
            .Line(Points<>{{getParameter(gamma_mm_lo),hist<>(ggmm.Transponate()).right().X().max()*1.5},{getParameter(gamma_mm_lo),0.0},
                           {getParameter(gamma_mm_hi),0.0},{getParameter(gamma_mm_hi),hist<>(ggmm.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '6γ missing mass, GeV/c^2'"<< "set xrange [2.2:3.4]"<<"set ylabel 'Events, n.d.'";
            const auto ggim=Hist(DATA, "All", histpath_central_reconstr, "GIM6");
            Plot("He36g-6gim-data",5)
            .Hist(ggim).Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM7"))
            .Line(Points<>{{getParameter(gamma_im_lo6),hist<>(ggim.Transponate()).right().X().max()*1.5},{getParameter(gamma_im_lo6),0.0},
                           {getParameter(gamma_im_hi6),0.0},{getParameter(gamma_im_hi6),hist<>(ggim.Transponate()).right().X().max()*1.5}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel '6γ invariant mass - Q_{3Heη}, GeV/c^2'"<< "set xrange [0.0:1.0]";
            const auto tim_hist=Hist(DATA, "All", histpath_central_reconstr, "TIM7-AllBins").Scale(4);
            Plot("He36g-tim-data",5)
                .Hist(tim_hist)
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xrange [-0.2:0.15]" << "set xlabel 'IM(^3He+6γ)-IM(p+d), GeV/c^2'";

            const auto DT=Hist(DATA, "All", histpath_central_reconstr, "dt2");
            const auto T=Hist(DATA, "All", histpath_central_reconstr, "t2");
            Plot("He36g-dt-data",5)
            .Hist(DT)
            .Line(Points<>{{getParameter(time_dt),0.0},{getParameter(time_dt),hist<>(DT.Transponate()).right().X().max()*1.5}},"condition")
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel 'dt γ-γ, ns'"<< "set key on";
            Plot("He36g-t-data",5)
            .Hist(T).Line(Points<>{{getParameter(time_t1),hist<>(T.Transponate()).right().X().max()*1.5},{getParameter(time_t1),0.0},
                {getParameter(time_t2),0.0},{getParameter(time_t2),hist<>(T.Transponate()).right().X().max()*1.5}},"condition")
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel 'dt ^3He-γ, ns'"<< "set key on";

            Plot("He36g-dt-data-final",5)
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt7"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel 'dt γ-γ, ns'"<< "set key on";
            Plot("He36g-t-data-final",5)
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t7"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"<<"set ylabel 'Events, n.d.'"
                    << "set xlabel 'dt ^3He-γ, ns'"<< "set key on";
    }
    ext_hist<2> ev_am,b_acc,ev_norm;
    vector<ext_hist<2>> acc;
    for (size_t i = 0; i < reaction.size(); i++) {
        acc.push_back(ext_hist<2>());
    }
    const auto lum_b_z = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc_z"));
    const auto lum_b_p = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc_p"));
    const auto lum_b_m = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc_m"));
    const list<size_t> params{pbeam_corr,he3_cut_h,he3_theta_cut,
        gamma_E_thr,time_dt,time_t1,time_t2,eta_theta_thr,he3mm_cut,
        gamma_mm_lo,gamma_mm_hi,gamma_im_lo6,gamma_im_hi6,three_pi0
    };
    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
                            << "Q є [" << setprecision(3)
                            << Q.min() << "; " << Q.max() << "] MeV").str();
        cout<<Qmsg << " plots3 & events count "<<endl;
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            cout<<Qmsg << " plots2 "<<r<<endl;
            acc[i]<<make_point(Q,RawSystematicError(params,[bin_num,&r,&histpath_central_reconstr](const string&suffix){
                const auto MC_TIM=Hist(MC, r, histpath_central_reconstr, "TIM7-Bin-"+to_string(bin_num),suffix).XRange(-0.3,0.3);
                const auto &N = NORM(r+suffix,[&histpath_central_reconstr,&r,&suffix](){
                    return Hist(MC, r, histpath_central_reconstr, "0-Reference",suffix);
                })[bin_num].Y();
                if (N.Above(0)) return extend_value<2,2>(std_error(MC_TIM.TotalSum().val())/N);
                else return uncertainties(0.,0.,0.);
            })());
        }
        cout<<Qmsg << " acceptance "<<endl;
        b_acc<<make_point(Q,SystematicError<bound_state_reaction_index>([&acc](const int i){return acc[i].right().Y();})());
        cout<<Qmsg << " events count "<<endl;
        ev_am<<make_point(Q,RawSystematicError(params,[bin_num,&histpath_central_reconstr](const string&suffix){
            const auto TIM=Hist(DATA, "All", histpath_central_reconstr, string("TIM7-Bin-") + to_string(bin_num),suffix).XRange(-0.3,0.3);
            return extend_value<1,2>(std_error(TIM.TotalSum().val()));
        })());
        cout<<Qmsg << " events count norm"<<endl;
        ev_norm<<make_point(Q,RawSystematicError(params,[
            bin_num,&histpath_central_reconstr,&reaction,
            &lum_b_m,&lum_b_p,&lum_b_z
        ](const string&suffix){
            const auto ac=SystematicError<bound_state_reaction_index>([&histpath_central_reconstr,&bin_num,&reaction,&suffix](const int i){
                const auto &r = reaction[i];
                const auto MC_TIM=Hist(MC, r, histpath_central_reconstr, "TIM7-Bin-"+to_string(bin_num),suffix).XRange(-0.3,0.3);
                static Cache<string,hist<>> NORM;
                const auto &N = NORM(r+suffix,[&histpath_central_reconstr,&r,&suffix](){
                    return Hist(MC, r, histpath_central_reconstr, "0-Reference",suffix);
                })[bin_num].Y();
                return extend_value<2,2>(std_error(MC_TIM.TotalSum().val())/N);
            })();
            const auto&lum=(suffix=="00+")?lum_b_p:(suffix=="00-")?lum_b_m:lum_b_z;
            const auto TIM=Hist(DATA, "All", histpath_central_reconstr, string("TIM7-Bin-") + to_string(bin_num),suffix).XRange(-0.3,0.3);
            return (extend_value<1,2>(std_error(TIM.TotalSum().val()))*trigger_he3_forward.scaling)/(ac*lum[bin_num].Y());
        })());
    }
    cout<<"Final plots"<<endl;
    Plot accplot("He36g-acceptance",5);
    accplot << "set title 'pd->^3He6γ'"<<"set key left top">>"set key right top"
            << "set xlabel 'Q_{3Heη}, MeV'"<<"set xrange [-70:30]"
            << "set ylabel 'Efficiency, n.d.'"
            << "set yrange [0:0.12]" << "set key on";
    accplot.Hist(wrap_hist(b_acc),rname[1]);
    for (size_t i = 3; i < reaction.size(); i++) {
        accplot.Hist(wrap_hist(acc[i]), rname[i]);
    }
    const auto he3eta_events = luminosity_he*branching_ratio*acc[3].XRange(10,30)
        *extend_hist<2,2>(hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed")).XRange(10,30))/trigger_he3_forward.scaling;

    Plot("He36g-events",5)
    .Hist(wrap_hist(ev_am),"data").Hist(wrap_hist(he3eta_events),"pd->^3Heη")
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:30]"
            << "set ylabel 'Events, n.d.'" << "set yrange [0:]"
            << "set title 'pd->^3He6γ " + runmsg + "'"<<"set key left top";
    const auto S2=take_uncertainty_component<1>(ev_am.XRange(10,30)-he3eta_events).TotalSum();
    Plot("He36g-events-BG",5)
    .Hist(hist<>(Points<value<>>{{S/1000.,2},{S2/1000.,1}}))
            <<"set yrange [0:3]"<<"set ytics ('right' 1,'left' 2)">>"unset ytics"
            << "set xlabel 'Events, 10^3'" << "set xrange [0:"+to_string(S2.val()*2./1000.)+"]"
            << "set title 'pd->^3He6γ background estimation " + runmsg + "'"<<"set key off";
    Plot("He36g-events-norm-bound",5)
        .Hist_2bars<1,2>(ev_norm.XRange(-70,10),"Statistical","Systematic","curve_3he_6gamma")
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
            << "set ylabel 'Normalized events, nb'" << "set yrange [0:60]"
            << "set title 'pd->^3He6γ "+runmsg+"'"<<"set key right top";
    Plot("He36g-events-norm-light",5)
        .Hist(wrap_hist(ev_norm).XRange(-70,10))
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q_{3Heη}, MeV'" << "set key on"<<"set xrange [-70:10]"
            << "set ylabel 'Normalized events, nb'" << "set yrange [0:60]"
            << "set title 'pd->^3He6γ "+runmsg+"'"<<"set key right top";
}
