// this file is distributed under 
// GPL v 3.0 license
#include <list>
#include <string>
#include <sstream>
#include <random>
#include <PBeamSmearing.h>
#include <PReaction.h>
#include <str_get.h>
#include <math_h/gnuplot/gnuplot.h>
#include "phys_constants.h"
using namespace std;
int main(int argc, char **arg){
	if(argc<3){
		printf("reaction expected\n");
		return -1;
	}
	printf("%s\n",arg[1]);
	string react="";
	for(int i=2;i<argc;i++){
		react+=string(arg[i]);
		if(i<(argc-1))react+=" ";
	}
	printf("%s\n",react.c_str());
	PUSHD();
	CD(PLUTO);
	PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
	smear->SetReaction("p+d");
	if(string(arg[1])=="all")
		smear->SetMomentumFunction(new TF1("Uniform","1",p_beam_low,p_beam_hi));
	if(string(arg[1])=="over")
		smear->SetMomentumFunction(new TF1("Uniform","1",p_he3_eta_threshold,p_beam_hi));
	makeDistributionManager()->Add(smear);
	std::default_random_engine gen;
	std::uniform_int_distribution<int> d(1,254);
	PUtils::SetSeed(d(gen));
	PReaction my_reaction(p_beam_hi,"p","d",
		const_cast<char*>(react.c_str()),
		const_cast<char*>(ReplaceAll(ReplaceAll(ReplaceAll(react," ",""),"[","_"),"]","_").c_str())
		,1,0,0,0);
	my_reaction.Loop(5000000);
	POPD();
	return 0;
}