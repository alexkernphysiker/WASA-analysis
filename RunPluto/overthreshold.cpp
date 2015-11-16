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
int main(int argc, char **arg){
#include "outpath.cc"
	if(argc<2){
		  printf("reaction expected\n");
		  return -1;
	}
	string react="";
	for(int i=1;i<argc;i++){
		react+=string(arg[i]);
		if(i<(argc-1))react+=" ";
	}
	printf("%s\n",react.c_str());
	string old=getcwd(NULL,0);
	chdir(outpath.c_str());
	PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
	smear->SetReaction("p+d");
	smear->SetMomentumFunction(new TF1("Uniform","1",p_he3_eta_threshold,p_beam_hi));
	makeDistributionManager()->Add(smear);
	std::default_random_engine gen;
	std::uniform_int_distribution<int> d(1,254);
	PUtils::SetSeed(d(gen));
	PReaction my_reaction(p_beam_hi,"p","d",
		const_cast<char*>(react.c_str()),
		const_cast<char*>(ReplaceAll(ReplaceAll(ReplaceAll(react," ",""),"[","_"),"]","_").c_str())
	,1,0,0,0);
	my_reaction.Loop(2000000);
	chdir(old.c_str());
	return 0;
}
