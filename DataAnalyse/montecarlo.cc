// this file is distributed under 
// GPL v 3.0 license
#include "montecarlo.h"
#include "../General/phys_constants.h"
MonteCarlo::MonteCarlo():Analysis(){
	AddLogSubprefix("MonteCarlo");
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();
	fMCVertexBank = MCTrf->GetVertexBank();
	fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));
}
MonteCarlo::~MonteCarlo(){}
double MonteCarlo::EventWeight(){
	return fEventHeader->GetWeight();
}
bool MonteCarlo::EventProcessingCondition(){
	return 
	gWasa->IsAnalysisMode(Wasa::kMCRaw)||
	gWasa->IsAnalysisMode(Wasa::kMCReco)||
	gWasa->IsAnalysisMode(Wasa::kMC);
}
double MonteCarlo::PBeam(){
	TVector3 result;
	result.SetMagThetaPhi(0,0,0);
	WVertexIter iterator(fMCVertexBank);
	int NrVertex=0;
	while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
		NrVertex++;
		if(NrVertex==2)
			for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
				WParticle *particle=vertex->GetParticle(particleindex);
				auto ptype=particle->GetType();
				for(auto P:first_particles)
					if(ptype==P.first){
						double ekin=particle->GetEkin();
						double theta=particle->GetTheta();
						double phi=particle->GetPhi();
						double m=P.second;
						double p=sqrt(ekin*(ekin+2*m));
						TVector3 P_vec;
						P_vec.SetMagThetaPhi(p,theta,phi);
						result=result+P_vec;
					}
			}
	}
	double res=result.Mag();
	if(res<=0)Log()<<"B_beam<=0";
	return res;
}
bool MonteCarlo::GetTrueParameters(ParticleType type,double&Ekin,double&theta,double&phi){
	WVertexIter iterator(fMCVertexBank);
	if(WVertex *vertex1=dynamic_cast<WVertex*>(iterator.Next()))
		if(WVertex *vertex2=dynamic_cast<WVertex*>(iterator.Next()))
			for(int particleindex=0,N=vertex2->NumberOfParticles();particleindex<N;particleindex++){
				WParticle *particle=vertex2->GetParticle(particleindex);
				auto ptype=particle->GetType();
				if(type==ptype){
					Ekin=particle->GetEkin();
					theta=particle->GetTheta();
					phi=particle->GetPhi();
					return true;
				}
			}
	Log(LogError)<<"cannot obtain true particle parameters";
	return false;
}
