// this file is distributed under 
// MIT license
#include <list>
#include <string>
#include <sstream>
#include <random>
#include <PBeamSmearing.h>
#include <PReaction.h>
#include <str_get.h>
#include <experiment_conv.h>
#include <reactions.h>
#include "config.h"
using namespace std;
string ReplaceAll(const string&source, const string& from, const string& to) {
	string str=source;
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return str;
}
string ReplaceAll(string&& str, const string& from, const string& to){
	return ReplaceAll(str,from,to);
}
int main(int argc, char **arg){
	if(argc<3){
		cout<<"reaction expected"<<endl;
		return -1;
	}
	string react="";
	for(int i=2;i<argc;i++){
		react+=string(arg[i]);
		if(i<(argc-1))react+=" ";
	}
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
	my_reaction.Loop(3000000);
	POPD();
	return 0;
}