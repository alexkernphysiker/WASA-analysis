// this file is distributed under 
// MIT license
#include "montecarlo.h"
MonteCarlo::MonteCarlo():Analysis(){
	AddLogSubprefix("MonteCarlo");
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();
	fMCVertexBank = MCTrf->GetVertexBank();
	fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));
}
MonteCarlo::~MonteCarlo(){}
void MonteCarlo::PrepairForEventAnalysis(){
	TVector3 result;
	result.SetMagThetaPhi(0,0,0);
	WVertexIter iterator(fMCVertexBank);
	if(dynamic_cast<WVertex*>(iterator.Next()))
		if(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next()))
			for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
				WParticle *particle=vertex->GetParticle(particleindex);
				ForFirstVertex([this,&result,particle](ParticleType type,double mass,Kinematic&out){
					if(type==particle->GetType()){
						out.E=particle->GetEkin();
						out.Th=particle->GetTheta();
						out.Phi=particle->GetPhi();
						double p=sqrt(out.E*(out.E+2*mass));
						TVector3 P_vec;
						P_vec.SetMagThetaPhi(p,out.Th,out.Phi);
						result=result+P_vec;
					}
				});
			}
	CachePBeam(result.Mag());
}
bool MonteCarlo::EventProcessingCondition(){
	return 
	gWasa->IsAnalysisMode(Wasa::kMCRaw)||
	gWasa->IsAnalysisMode(Wasa::kMCReco)||
	gWasa->IsAnalysisMode(Wasa::kMC);
}
