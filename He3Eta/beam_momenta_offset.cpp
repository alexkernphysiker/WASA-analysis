// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/functions.h>
#include <math_h/sigma.h>
#include <math_h/interpolate.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Kinematics/particles.h>
#include <Kinematics/reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"beam_momenta");
	Reaction He3eta(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	LinearInterpolation<double> Q2P=SortedPoints<double>(
		[&He3eta](double p){return He3eta.P2Q(p);},ChainWithStep(0.0,0.0002,3.0)
	).Transponate();
	LinearInterpolation<double> Q2E=SortedPoints<double>(
		[&He3eta](double e){return He3eta.E2Q(e);},ChainWithStep(0.0,0.0002,3.0)
	).Transponate();
	LinearInterpolation<double> ThetaMax2P=SortedPoints<double>(
		[&He3eta](double p){
			const double Ek_theta_max=0.3;
			return He3eta.PbEr2Theta(p,Ek_theta_max)*180/PI();
			
		},ChainWithStep(He3eta.PThreshold(),0.001,3.0)
	).Transponate();
	string he3eta="He3eta";
	SortedPoints<value<double>> offs_mc,offs_data,ek_mc,ek_data,ek_before_shift;
	auto QBins=Hist(MC,he3eta,{"Histograms","He3Forward_Reconstruction"},"0-Reference");
	for(size_t bin_num=9,bin_count=QBins.size()-1;bin_num<bin_count;bin_num++){
		auto Q=QBins[bin_num].X();
		double p=Q2P(Q.val()/1000.0);
		auto P_offset=[&He3eta,&ThetaMax2P,&Q,&p](const hist2d<double>&kin_,const double ratio){
			double max=0;StandardDeviation<double> theta_avr;
			vector<point<value<double>>> points;
			kin_.FullCycle([&max](const point3d<value<double>>&P){
				if(max<P.Z().val())max=P.Z().val();
			});
			kin_.FullCycle([&max,&theta_avr,&points,&ratio](const point3d<value<double>>&P){
				if(
					(P.X().val()>0.2)&&(P.X().val()<0.4)
					&&(P.Y().val()<7.3)
					&&(P.Z().val()>(ratio*max))
				){
					points.push_back(point<value<double>>(P.X(),P.Y()));
					for(unsigned long i=0;i<(P.Z().val());i++)
						theta_avr<<P.Y().val();
				}
			});
			auto P=func_value(ThetaMax2P.func(),theta_avr());
			Plot<double>().Hist(points)
			.Line(SortedPoints<double>(
				[&He3eta,&theta_avr,&P](double E)->double{return He3eta.PbEr2Theta(P.val(),E)*180./PI();},
				ChainWithStep(0.25,0.001,0.35))
			)
			<<"set xrange [0.2:0.4]"<<"set yrange [4.5:7.5]";
			return point<value<double>>(Q,P-value<double>(p));
		};
		auto Ek_avr=[&He3eta,&ThetaMax2P,&Q,&p](const hist2d<double>&kin_,const double ratio){
			double max=0;StandardDeviation<double> ek_avr;
			vector<point<value<double>>> points;
			kin_.FullCycle([&max](const point3d<value<double>>&P){
				if(max<P.Z().val())max=P.Z().val();
			});
			kin_.FullCycle([&max,&ek_avr,&points,&ratio](const point3d<value<double>>&P){
				if(
					(P.X().val()>0.2)&&(P.X().val()<0.4)
					&&(P.Y().val()<7.3)
					&&(P.Z().val()>(ratio*max))
				){
					points.push_back(point<value<double>>(P.X(),P.Y()));
					for(unsigned long i=0;i<(P.Z().val()-(max*ratio));i++)
						ek_avr<<P.X().val();
				}
			});
			return point<value<double>>(Q,ek_avr());
		};
		auto kin_v=Hist2d(MC,he3eta,{"Histograms","He3Forward_Vertices"},string("Kinematic-vertex-Bin-")+to_string(bin_num)).Scale(3,3);
		PlotHist2d<double>(sp2).Distr(kin_v) << "set xlabel 'E_k, GeV'" << "set ylabel 'theta, deg'"
		<< "set xrange [0.2:0.4]";
		auto kin_mc=Hist2d(MC,he3eta,{"Histograms","He3Forward_Reconstruction"},string("Kinematic-reconstructed-Bin-")+to_string(bin_num)).Scale(3,3);
		PlotHist2d<double>(sp2).Distr(kin_mc)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'"
		<< "set xrange [0.2:0.4]";
		auto data_hist_init=Hist2d(DATA,"He3",{"Histograms","He3Forward_Reconstruction"},string("Kinematic-reconstructed-Bin-")+to_string(bin_num));
		{
			hist2d<double> shifted_back(
				BinsByCount(data_hist_init.X().size(),
					    data_hist_init.X().left().min()-he3_forward_correct_energy,
					    data_hist_init.X().right().max()-he3_forward_correct_energy
				),
				BinsByCount(data_hist_init.Y().size(),
					    data_hist_init.Y().left().min(),
					    data_hist_init.Y().right().max()
				)
			);
			for(size_t x=0;x<shifted_back.X().size();x++)
				for(size_t y=0;y<shifted_back.Y().size();y++)
					shifted_back.Bin(x,y)=data_hist_init[x][y];
			shifted_back=shifted_back.Scale(3,3);
			PlotHist2d<double>(sp2).Distr(shifted_back)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'"
			<< "set xrange [0.2:0.4]";
			ek_before_shift << Ek_avr(shifted_back,3.0/4.0);
		}
		auto kin_data=data_hist_init.Scale(3,3);
		PlotHist2d<double>(sp2).Distr(kin_data)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'"
		<< "set xrange [0.2:0.4]";
		offs_mc << P_offset(kin_mc,3.0/4.0);
		offs_data << P_offset(kin_data,3.0/4.0);
		ek_mc << Ek_avr(kin_mc,3.0/4.0);
		ek_data << Ek_avr(kin_data,3.0/4.0);
	}
	Plot<double>().Hist(offs_mc,"WMC").Hist(offs_data,"Data")
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'delta P, GeV/c'";
	Plot<double>().Hist(ek_mc,"WMC").Hist(ek_data,"Data")
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'Ek_{avr}, GeV'";
	return 0;
}