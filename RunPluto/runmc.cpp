#include <list>
#include <string>
#include <PBeamSmearing.h>
#include <PReaction.h>
#include "math_h/randomfunc.h"
const double he3_eta_threshold=1.5727;
const double beam_low=1.426;
const double beam_hi=1.635;
using namespace std;
string ReplaceAll(string str, const string& from, const string& to) {
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return str;
}
int main(int, char **){
	PUtils::SetSeed(RandomUniformlyI(1,50));
	PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
	smear->SetReaction("p+d");
	smear->SetMomentumFunction(new TF1("Uniform","1",he3_eta_threshold,beam_hi));
	makeDistributionManager()->Add(smear);
	list<string> reactlist;
	reactlist.push_back("He3 eta");
	reactlist.push_back("He3 eta [g g]");
	reactlist.push_back("He3 eta [pi0 pi0 pi0]");
	reactlist.push_back("p p p pi-");
	reactlist.push_back("d p pi0");
	reactlist.push_back("p p n pi0");
	for(auto react:reactlist){
		PReaction my_reaction(beam_hi,"p","d",
			const_cast<char*>(react.c_str()),
			const_cast<char*>(ReplaceAll(ReplaceAll(ReplaceAll(react," ",""),"[","_"),"]","_").c_str())
		,1,0,0,0);
		my_reaction.Loop(1000000);
	}
	return 0;
}
