#include <list>
#include <string>
#include <PBeamSmearing.h>
#include <PReaction.h>
#include <phys_constants.h>
#include "math_h/randomfunc.h"
#include <replace.cc>
using namespace std;
int main(int, char **){
	PUtils::SetSeed(RandomUniformlyI(1,50));
	PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
	smear->SetReaction("p+d");
	smear->SetMomentumFunction(new TF1("Uniform","1",p_he3_eta_threshold,p_beam_hi));
	makeDistributionManager()->Add(smear);
	list<string> reactlist;
	reactlist.push_back("He3 eta");
	reactlist.push_back("He3 eta [g g]");
	reactlist.push_back("He3 eta [pi0 pi0 pi0]");
	reactlist.push_back("p p p pi-");
	reactlist.push_back("d p pi0");
	reactlist.push_back("p p n pi0");
	for(auto react:reactlist){
		PReaction my_reaction(p_beam_hi,"p","d",
			const_cast<char*>(react.c_str()),
			const_cast<char*>(ReplaceAll(ReplaceAll(ReplaceAll(react," ",""),"[","_"),"]","_").c_str())
		,1,0,0,0);
		my_reaction.Loop(1000000);
	}
	return 0;
}
