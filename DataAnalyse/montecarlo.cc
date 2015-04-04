#include "montecarlo.h"
MonteCarlo::MonteCarlo():Analysis(){
	WTrackFinder *MCTrf = dynamic_cast<WTrackFinder*>(gDataManager->GetAnalysisModule("MCTrackFinder","default"));
	fMCTrackBank  = MCTrf->GetTrackBank();
	fMCVertexBank = MCTrf->GetVertexBank();
	fEventHeader = dynamic_cast<REventWmcHeader*>(gDataManager->GetDataObject("REventWmcHeader","EventHeader"));
}
MonteCarlo::~MonteCarlo(){}
double MonteCarlo::EventWeight(){return fEventHeader->GetWeight();}
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
	return result.Mag();
}
MonteCarlo::CheckHists::CheckHists(ParticleType t){
	type=t;
	Ekin=new TH1F(Form("Check_E_%i",int(t)),"",500,-0.1,0.1);
	Theta=new TH1F(Form("Check_Theta_%i",int(t)),"",500,-1,1);
	Phi=new TH1F(Form("Check_Phi_%i",int(t)),"",500,-1,1);
}
void MonteCarlo::PrepareCheck(){
	Analysis::PrepareCheck();
	for(auto P:first_particles){
		CheckHists h(P.first);
		check.push_back(h);
		gHistoManager->Add(h.Ekin,"Reconstruction");
		gHistoManager->Add(h.Theta,"Reconstruction");
		gHistoManager->Add(h.Phi,"Reconstruction");
	}
}
void MonteCarlo::CheckParticleTrack(ParticleType type, double Ekin, double theta, double phi){
	const double two_pi=2*3.1415926;
	while(phi<0)phi+=two_pi;
	while(phi>=two_pi)phi-=two_pi;
	WVertexIter iterator(fMCVertexBank);
	int NrVertex=0;
	for(CheckHists H:check)
		if(H.type==type)
			while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
				NrVertex++;
				if(NrVertex==2)
					for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
						WParticle *particle=vertex->GetParticle(particleindex);
						auto ptype=particle->GetType();
						if(type==ptype){
							double p_ekin=particle->GetEkin();
							double p_theta=particle->GetTheta();
							double p_phi=particle->GetPhi();
							while(p_phi<0)p_phi+=two_pi;
							while(p_phi>=two_pi)p_phi-=two_pi;
							H.Ekin->Fill(Ekin-p_ekin);
							H.Theta->Fill(theta-p_theta);
							H.Phi->Fill(phi-p_phi);
						}
					}
			}
}
