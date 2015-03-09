#include <list>
#include <string>
#include <PBeamSmearing.h>
#include <PReaction.h>
using namespace std;
int main(int, char **){
	PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam smearing");
	smear->SetReaction("p + d");
	smear->SetMomentumFunction(new TF1("Uniform","1",1.426,1.635));
	makeDistributionManager()->Add(smear);
	list<string> reactlist;
	reactlist.push_back("He3 eta [g g]");
	reactlist.push_back("He3 eta [pi0 pi0 pi0]");
	for(string react:reactlist){
		PReaction my_reaction(1.635,"p","d",
			const_cast<char*>(react.c_str()), const_cast<char*>((react).c_str())
		,1,0,0,0);
		my_reaction.Loop(1000000);
	}
	return 0;
}
