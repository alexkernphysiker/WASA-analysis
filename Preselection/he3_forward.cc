// this file is distributed under 
// MIT license
#include <TLorentzVector.h>
#include <CAnalysisModule.hh>
#include <CDataManager.hh>
#include <FDEdep2Ekin.hh>
#include <CCardWDET.hh>
#include <Wasa.hh>
#include <CAnalysisModule.hh>
#include <REventWmcHeader.hh>
#include <REventHeader.hh>
#include <WTrackBank.hh>
#include <WVertexBank.hh>
#include <FDFTHTracks.hh>
#include <CDTracksSimple.hh>
#include <TCutG.h>
#include "../phys_constants.h"
#include "../kinematics.h"
#include "../reconstruction_types.h"
#include "trackprocessing.h"
#include "detectors.h"
#include "reconstruction.h"
#include "reactions.h"
#include "data.h"
#include "montecarlo.h"
namespace ReactionSetup{
	using namespace std;
	using namespace TrackAnalyse;
	string dirname(){return "He3Forward_Reconstruction";};
	Analysis* Prepare(He3Modification mode){
		Analysis* res=nullptr;
		switch(mode){
			case forData:
				res=new RealData();
				break;
			case forEta:
				res=new MonteCarlo();
				res->AddParticleToFirstVertex(kHe3,m_3He);
				res->AddParticleToFirstVertex(kEta,m_eta);
				break;
			case forPi0:
				res=new MonteCarlo();
				res->AddParticleToFirstVertex(kHe3,m_3He);
				res->AddParticleToFirstVertex(kPi0,m_pi0);
				break;
		};
		return res;
	}
	shared_ptr<AbstractChain> ReconstructionProcess(const Analysis&data,const Axis&Q){
		return make_shared<ChainCheck>()
		<<Forward::Get().CreateMarker(dirname(),"1-AllTracks")
		<<make_shared<Hist1D>(dirname(),"1-AllTracks",Q)
		<<[](WTrack&T)->bool{
			//ToDo: replace by more reasonable condition
			return (T.Theta()!=0.125);
		}
		<<Forward::Get().CreateMarker(dirname(),"2-FPC")<<make_shared<Hist1D>(dirname(),"2-FPC",Q)
		<<(make_shared<ChainOr>()//E_dep cuts
			<<(make_shared<ChainCheck>()
				<<[](WTrack&T)->bool{return Forward::Get().StoppingLayer(T)==kFRH1;}
				<<[](WTrack&T)->bool{
					//Achtung - static
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("FRH1_cut",15);
						cut->SetVarX("FRH1");
						cut->SetVarY("FTH1");
						cut->SetPoint(15,0.019,0.025);
						cut->SetPoint(14,0.069,0.019);
						cut->SetPoint(13,0.121,0.016);
						cut->SetPoint(12,0.162,0.014);
						cut->SetPoint(11,0.206,0.013);
						cut->SetPoint(10,0.304,0.012);
						cut->SetPoint(9,0.4,0.012);
						cut->SetPoint(8,0.4,0.008);
						cut->SetPoint(7,0.298,0.008);
						cut->SetPoint(6,0.201,0.010);
						cut->SetPoint(5,0.141,0.012);
						cut->SetPoint(4,0.105,0.013);
						cut->SetPoint(3,0.061,0.017);
						cut->SetPoint(2,0.027,0.021);
						cut->SetPoint(1,0.019,0.025);
					}
					double x=Forward::Get()[kFRH1].Edep(T);
					double y=Forward::Get()[kFTH1].Edep(T);
					return cut->IsInside(x,y);
				}
				<<Forward::Get().CreateMarker(dirname(),"2.5-FRH1")
				<<make_shared<Parameter>([&data](WTrack&track)->double{
					//Achtung - static
					static FitBasedReconstruction<Reconstruction::He3EnergyFRH1,WTrack&> energy(
						"He3.E.FRH1",
						{
							[](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},
							[](WTrack&track){return track.Theta();}
						},
						[&data](WTrack&){return data.FromFirstVertex(kHe3).E;}
					);
					return energy.Reconstruct(track);
				})
			)
			<<(make_shared<ChainCheck>()
				<<[](WTrack&T)->bool{return Forward::Get().StoppingLayer(T)==kFRH2;}
				<<[](WTrack&T)->bool{
					vector<double> Ed={
						Forward::Get()[kFRH1].Edep(T),
						Forward::Get()[kFRH2].Edep(T)
					};
					return (Ed[0]>(0.25-0.417*Ed[1]))&&(Ed[0]<(0.35-0.417*Ed[1]))&&(Ed[1]<0.22);
				}
				<<Forward::Get().CreateMarker(dirname(),"2.5-FRH2")
				<<make_shared<Parameter>([&data](WTrack&track)->double{
					//Achtung - static
					static FitBasedReconstruction<Reconstruction::He3EnergyFRH2,WTrack&> energy(
						"He3.E.FRH2",
						{
							[](WTrack&T){return Forward::Get()[kFRH1].Edep(T)+Forward::Get()[kFRH2].Edep(T);},
														    [](WTrack&T){return T.Theta();}
						},
						[&data](WTrack&){return data.FromFirstVertex(kHe3).E;}
					);
					return energy.Reconstruct(track);
				})	
			)
		)//end E_dep cuts
		<<Forward::Get().CreateMarker(dirname(),"3-AllCuts")<<make_shared<Hist1D>(dirname(),"3-AllCuts",Q)
		<<make_shared<Parameter>([](WTrack&T)->double{return T.Theta();})
		<<make_shared<Parameter>([](WTrack&T)->double{return T.Phi();})
		<<[](const vector<double>&P)->bool{return isfinite(P[0])&&isfinite(P[1])&&isfinite(P[2]);}
		<<Forward::Get().CreateMarker(dirname(),"4-Reconstructed")
		<<make_shared<Hist1D>(dirname(),"4-Reconstructed",Q);
	}
	shared_ptr<AbstractChain> MissingMass(const Analysis&data,const Axis&Q){
		Axis M([](const vector<double>&P)->double{return P[3];},m_pi0-0.5,m_eta+0.5,200);
		return make_shared<Chain>()
		<<make_shared<Parameter>([&data](const vector<double>&He3)->double{
			double p=sqrt(He3[0]*(He3[0]+2*m_3He));
			TVector3 p_He3;
			p_He3.SetMagThetaPhi(p,He3[1],He3[2]);
			TLorentzVector P_He3;
			P_He3.SetVectM(p_He3,m_3He);
			TLorentzVector P_Total;{
				TVector3 p_beam;
				p_beam.SetMagThetaPhi(data.PBeam(),0,0);
				TLorentzVector P_Beam;
				TLorentzVector P_Target;
				P_Beam.SetVectM(p_beam,m_p);
				TVector3 ptarget;
				ptarget.SetMagThetaPhi(0,0,0);
				P_Target.SetVectM(ptarget,m_d);
				P_Total=P_Beam+P_Target;
			}
			TLorentzVector P_Missing=P_Total-P_He3;
			return P_Missing.M();
		})
		<<make_shared<SetOfHists1D>(dirname(),"MissingMass",Q,M);
	}
	Axis Ek([](const vector<double>&P)->double{return P[0];},0.1,0.6,100);
	Axis Th([](const vector<double>&P)->double{return P[1];},0.06,0.16,100);
	Axis Phi([](const vector<double>&P)->double{return P[2];},0.00,6.28,360);
	shared_ptr<AbstractChain> He3Eta_cut(const Analysis&data){
		//The axis is duplicated because cuts are strongly connected with it
		return make_shared<ChainBinner>(Axis([&data]()->double{return 1000.0*Q_He3eta(data.PBeam());},0.0,30.0,12))
				<<[]()->bool{return false;}//0
				<<[]()->bool{return false;}
				<<[]()->bool{return false;}//2
				<<[]()->bool{return false;}
				<<[]()->bool{return false;}//4
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut5",5);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1,0.279,0.095);
						cut->SetPoint(2,0.279,0.091);
						cut->SetPoint(3,0.323,0.091);
						cut->SetPoint(4,0.323,0.095);
						cut->SetPoint(5,0.279,0.095);
					}
					return cut->IsInside(Ek(T,P),Th(T,P));
				}
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut6",5);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1,0.264,0.091);
						cut->SetPoint(2,0.336,0.091);
						cut->SetPoint(3,0.310,0.104);
						cut->SetPoint(4,0.290,0.104);
						cut->SetPoint(5,0.264,0.091);
					}
					return cut->IsInside(Ek(T,P),Th(T,P));
				}//6
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut7",7);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1,0.263,0.091);
						cut->SetPoint(2,0.342,0.091);
						cut->SetPoint(3,0.332,0.104);
						cut->SetPoint(4,0.311,0.112);
						cut->SetPoint(5,0.286,0.112);
						cut->SetPoint(6,0.273,0.104);
						cut->SetPoint(7,0.263,0.091);
					}
					return cut->IsInside(Ek(T,P),Th(T,P));
				}
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut8",10);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1 ,0.259,0.091);
						cut->SetPoint(2 ,0.282,0.091);
						cut->SetPoint(3 ,0.303,0.098);
						cut->SetPoint(4 ,0.322,0.091);
						cut->SetPoint(5 ,0.354,0.091);
						cut->SetPoint(6 ,0.338,0.107);
						cut->SetPoint(7 ,0.317,0.115);
						cut->SetPoint(8 ,0.291,0.115);
						cut->SetPoint(9 ,0.267,0.107);
						cut->SetPoint(10,0.259,0.091);
					}
					return cut->IsInside(Ek(T,P),Th(T,P));
				}//8
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut9",10);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1 ,0.252,0.091);
						cut->SetPoint(2 ,0.278,0.091);
						cut->SetPoint(3 ,0.302,0.101);
						cut->SetPoint(4 ,0.329,0.091);
						cut->SetPoint(5 ,0.362,0.091);
						cut->SetPoint(6 ,0.349,0.110);
						cut->SetPoint(7 ,0.320,0.120);
						cut->SetPoint(8 ,0.286,0.120);
						cut->SetPoint(9 ,0.263,0.110);
						cut->SetPoint(10,0.252,0.091);
					}
					return cut->IsInside(Ek(T,P),Th(T,P));
				}
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut10",10);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1 ,0.249,0.091);
						cut->SetPoint(2 ,0.270,0.091);
						cut->SetPoint(3 ,0.304,0.108);
						cut->SetPoint(4 ,0.339,0.091);
						cut->SetPoint(5 ,0.368,0.091);
						cut->SetPoint(6 ,0.353,0.114);
						cut->SetPoint(7 ,0.318,0.125);
						cut->SetPoint(8 ,0.293,0.125);
						cut->SetPoint(9 ,0.263,0.114);
						cut->SetPoint(10,0.249,0.091);
					}
					return cut->IsInside(Ek(T,P),Th(T,P));
				}//10
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut11",11);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1 ,0.248,0.091);
						cut->SetPoint(2 ,0.268,0.091);
						cut->SetPoint(3 ,0.304,0.110);
						cut->SetPoint(4 ,0.343,0.091);
						cut->SetPoint(5 ,0.379,0.091);
						cut->SetPoint(6 ,0.365,0.114);
						cut->SetPoint(7 ,0.320,0.130);
						cut->SetPoint(8 ,0.292,0.130);
						cut->SetPoint(9 ,0.253,0.114);
						cut->SetPoint(10,0.248,0.091);
					}
					return cut->IsInside(Ek(T,P),Th(T,P));
				}
			;
	}
	shared_ptr<AbstractChain> KinematicHe3Test(const Analysis&data,const Axis&Q,bool MC){
		auto res=make_shared<Chain>()<<make_shared<SetOfHists2D>(dirname(),"Kinematic-reconstructed",Q,Ek,Th);
		if(MC){
			Axis Ev([&data]()->double{return data.FromFirstVertex(kHe3).E;},Ek);
			Axis Tv([&data]()->double{return data.FromFirstVertex(kHe3).Th;},Th);
			res<<make_shared<SetOfHists2D>(dirname(),"Kinematic-vertex",Q,Ev,Tv);
		}
		return res;
	}
	///Reaction analysis types visible from reactions.h
	Analysis* He3_forward_analyse(He3Modification mode){
		auto res=Prepare(mode);
		Axis Q([res]()->double{return 1000.0*Q_He3eta(res->PBeam());},0.0,30.0,12);
		if(forData!=mode)
			res->EventProcessing()<<make_shared<Hist1D>(dirname(),"0-Reference",Q);
		res->TrackTypeProcess(kFDC)<<(make_shared<ChainCheck>()
			<<ReconstructionProcess(*res,Q)
			<<He3Eta_cut(*res)<<MissingMass(*res,Q)
			<<KinematicHe3Test(*res,Q,mode==forEta)
		);
		return res;
	}
	Analysis* He3_forward_reconstruction(He3Modification mode){
		auto res=Prepare(mode);
		Axis Q([res]()->double{return 1000.0*Q_He3eta(res->PBeam());},0.0,30.0,12);
		if(forData!=mode)
			res->EventProcessing()<<make_shared<Hist1D>(dirname(),"0-Reference",Q);
		res->TrackTypeProcess(kFDC)<<(make_shared<ChainCheck>()
			<<ReconstructionProcess(*res,Q)
			<<KinematicHe3Test(*res,Q,mode)
		);
		return res;
	}
}
