// this file is distributed under 
// GPL v 3.0 license
#include <list>
#include <string>
#include <sstream>
#include <random>
#include <PBeamSmearing.h>
#include <PReaction.h>
#include <phys_constants.h>
#include "math_h/gnuplot/gnuplot.h"
using namespace std;
int main(int, char **){
	string outpath;{
		stringbuf buffer;
		ostream os (&buffer); 
		os<<getenv("PLUTO_OUTPUT");
		outpath=buffer.str();
		printf("output path: %s\n",outpath.c_str());
	}
	string old=getcwd(NULL,0);
	chdir(outpath.c_str());
	std::default_random_engine gen;
	std::uniform_int_distribution<int> d(1,254);
	list<pair<string,double>> reactlist;
	reactlist.push_back(make_pair("He3 eta",p_he3_eta_threshold));
	reactlist.push_back(make_pair("He3 pi0 pi0",p_beam_low));
	reactlist.push_back(make_pair("He3 pi0 pi0 pi0",p_beam_low));
	for(auto react:reactlist){
		PUtils::SetSeed(d(gen));
		PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
		smear->SetReaction("p+d");
		smear->SetMomentumFunction(new TF1("Uniform","1",react.second,p_beam_hi));
		makeDistributionManager()->Add(smear);
		PReaction my_reaction(p_beam_hi,"p","d",
			const_cast<char*>(react.first.c_str()),
			const_cast<char*>(ReplaceAll(ReplaceAll(ReplaceAll(react.first," ",""),"[","_"),"]","_").c_str())
		,1,0,0,0);
		my_reaction.Loop(1000000);
		delete smear;
	}
	chdir(old.c_str());
	return 0;
}
