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
int main(){
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_central_2gamma");
    vector<string> histpath_central_reconstr={"Histograms","CentralGammas2"};
    vector<string> reaction={"bound1-2g","He3eta","He3pi0","He3pi0pi0","He3pi0pi0pi0"};
    hist<double> norm=Hist(MC,reaction[0],{"Histograms","CentralGammas"},"0-Reference");
    hist<double> norm2=Hist(MC,reaction[1],{"Histograms","CentralGammas"},"0-Reference");
    const auto runs=PresentRuns("C");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";
    Plot<> theory,experiment;
    for(const auto&r:reaction){
	theory.Line(Hist(MC,r,histpath_central_reconstr,"GIMDiff-AllBins").toLine(),r);
    }
    theory<<"set key on"<< "set yrange [0:]";
    experiment.Hist(Hist(DATA,"C",histpath_central_reconstr,"GIMDiff-AllBins"),"DATA")
    <<"set key on"<<"set title '"+runmsg+"'"<< "set yrange [0:]";
    hist<double> ev_am,acceptance,acceptance2;
    for(size_t bin_num=0,bin_count=norm.size();bin_num<bin_count;bin_num++){
	const auto&Q=norm[bin_num].X();
	const auto&N=norm[bin_num].Y();
	const auto&N2=norm2[bin_num].Y();
	const string Qmsg=static_cast<stringstream&>(stringstream()
	    <<"Q in ["<<setprecision(3)
	    <<Q.min()<<"; "<<Q.max()<<"] MeV"
	).str();
	const hist<double> mc=Hist(MC,reaction[0],histpath_central_reconstr,string("GIMDiff-Bin-")+to_string(bin_num)),
	mc2=Hist(MC,reaction[1],histpath_central_reconstr,string("GIMDiff-Bin-")+to_string(bin_num)),
	data=Hist(DATA,"C",histpath_central_reconstr,string("GIMDiff-Bin-")+to_string(bin_num)),
	MC=mc.XRange(0,0.1),MC2=mc2.XRange(0,0.1),DATA=data.XRange(0,0.1);
	Plot<>().Hist(data).Hist(DATA)<<"set title '"+Qmsg+";"+runmsg+"'"<< "set yrange [0:]";
	Plot<>().Hist(mc).Hist(MC)<<"set title '"+Qmsg+";MC "+reaction[0]+"'"<< "set yrange [0:]";
	acceptance<<point<value<double>>(Q,MC.TotalSum()/N);
	if(Q>0.0){
	    Plot<>().Hist(mc2).Hist(MC2)<<"set title '"+Qmsg+";MC "+reaction[1]+"'"<< "set yrange [0:]";
	    acceptance2<<point<value<double>>(Q,MC2.TotalSum()/N2);
	}
	ev_am<<point<value<double>>(Q,DATA.TotalSum());
    }
    Plot<>().Hist(acceptance,reaction[0]).Hist(acceptance2,reaction[1])
    <<"set title 'Acceptance'"<< "set yrange [0:]"<<"set key on";
    Plot<>().Hist(ev_am)<<"set title 'Data events'"<< "set yrange [0:]";
    Plot<>().Hist((ev_am/acceptance)*trigger_he3_forward.scaling)<<"set title 'Events norm'"<< "set yrange [0:]";
}


