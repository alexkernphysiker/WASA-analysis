#include <memory>
#include "routines.h"
using namespace std;
TH1F *He3_Ekin, *He3_Theta, *He3_Phi;
void InitAnalysis(CHistoManager *histoManager){
	He3_Ekin=new TH1F("Kinetic Energy","",1000,0,2);
	histoManager->Add(He3_Ekin,"E_test");
	He3_Theta=new TH1F("Theta","",18,0.0,30.0);
	histoManager->Add(He3_Theta,"Theta_test");
	He3_Phi=new TH1F("Phi","",36, 0,360);
	histoManager->Add(He3_Phi,"Phi_test");
}
void ProcessVertexBank(WVertexBank *vertexBank){
	WVertexIter iterator(vertexBank);
	while(WVertex *vertex=dynamic_cast<WVertex*>(iterator.Next())){
		for(int particleindex=0; particleindex<vertex->NumberOfParticles(); particleindex++){
			WParticle *particle=vertex->GetParticle(particleindex);
			if(kHe3==particle->GetType()){
				He3_Ekin->Fill(particle->GetEkin());
				He3_Theta->Fill(particle->GetTheta()*180/3.1415926);
				He3_Phi->Fill(particle->GetPhi()*180/3.1415926);
			}
		}
	}
}
