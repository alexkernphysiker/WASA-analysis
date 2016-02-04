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
	Axis Q_axis(Analysis*res){return Axis([res]()->double{return 1000.0*Q_He3eta(res->PBeam());},0.0,30.0,12);}
	Axis Th_he([](const vector<double>&P)->double{return P[0]*180.0/3.1415926;},3.5,9.0,550);
	Axis Phi_he([](const vector<double>&P)->double{return P[1]*180.0/3.1415926;},0.0,360.0,360);
	Axis Ek_he([](const vector<double>&P)->double{return P[2];},0.1,0.6,500);
	Axis MM_he([](const vector<double>&P)->double{return P[3];},m_eta-0.02,m_eta+0.02,40);
	shared_ptr<AbstractChain> ReconstructionProcess(const Analysis&data,const Axis&Q){
		return make_shared<ChainCheck>()
		<<Forward::Get().CreateMarker(dirname(),"1-AllTracks")
		<<make_shared<Hist1D>(dirname(),"1-AllTracks",Q)
		<<[](WTrack&T)->bool{
			//ToDo: replace by more reasonable condition
			return (T.Theta()!=0.125);
		}
		<<make_shared<Parameter>([](WTrack&T)->double{return T.Theta();})
		<<make_shared<Parameter>([](WTrack&T)->double{return T.Phi();})
		<<Forward::Get().CreateMarker(dirname(),"2-FPC")<<make_shared<Hist1D>(dirname(),"2-FPC",Q)
		<<[](WTrack&T,const vector<double>&P)->bool{
			return (Th_he(T,P)<7.5);// Due to Kinematical predictions
		}
		<<(make_shared<ChainOr>()//E_dep cuts
			<<(make_shared<ChainCheck>()//for particles stopped in FRH1
				<<[](WTrack&T)->bool{return Forward::Get().StoppingLayer(T)==kFRH1;}
				<<[](WTrack&T)->bool{
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
				<<make_shared<Parameter>([&data](WTrack&track)->double{//Reconstructing kinetic energy
					static FitBasedReconstruction<Reconstruction::He3EnergyFRH1,WTrack&> energy(
						"He3.E.FRH1",{[](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},[](WTrack&track){return track.Theta();}},
						[&data](WTrack&){return data.FromFirstVertex(kHe3).E;}
					);
					return energy.Reconstruct(track);
				})
			)
			<<(make_shared<ChainCheck>()//for particles stopped in FRH2
				<<[](WTrack&T)->bool{return Forward::Get().StoppingLayer(T)==kFRH2;}
				<<[](WTrack&T)->bool{
					if(Forward::Get()[kFRH2].Edep(T)>0.22)return false;
					double locusline=0.3-0.417*Forward::Get()[kFRH2].Edep(T);
					return IsIn(Forward::Get()[kFRH1].Edep(T),make_pair(locusline-0.05,locusline+0.05));
				}
				<<Forward::Get().CreateMarker(dirname(),"2.5-FRH2")
				<<make_shared<Parameter>([&data](WTrack&track)->double{//Reconstructing kinetic energy
					static FitBasedReconstruction<Reconstruction::He3EnergyFRH2,WTrack&> energy(
						"He3.E.FRH2",{[](WTrack&T){return Forward::Get()[kFRH1].Edep(T)+Forward::Get()[kFRH2].Edep(T);},[](WTrack&T){return T.Theta();}},
						[&data](WTrack&){return data.FromFirstVertex(kHe3).E;}
					);
					return energy.Reconstruct(track);
				})	
			)
		)//end E_dep cuts
		<<Forward::Get().CreateMarker(dirname(),"3-AllCuts")<<make_shared<Hist1D>(dirname(),"3-AllCuts",Q)
		<<[](const vector<double>&P)->bool{return isfinite(P[0])&&isfinite(P[1])&&isfinite(P[2]);}
		<<Forward::Get().CreateMarker(dirname(),"4-Reconstructed")
		<<make_shared<Hist1D>(dirname(),"4-Reconstructed",Q);
	}
	shared_ptr<AbstractChain> MissingMass(const Analysis&data,const Axis&Q){
		return make_shared<Chain>()
		<<make_shared<Parameter>([&data](WTrack&T,const vector<double>&P)->double{
			double p=sqrt(Ek_he(T,P)*(Ek_he(T,P)+2*m_3He));
			TVector3 p_He3;
			p_He3.SetMagThetaPhi(p,Th_he(T,P),Phi_he(T,P));
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
		<<make_shared<SetOfHists1D>(dirname(),"MissingMass",Q,MM_he);
	}
	shared_ptr<AbstractChain> He3Eta_kin_cut(const Analysis&data,const Axis&Q){
		return make_shared<ChainCheck>()
			<<(make_shared<ChainBinner>(Q)
				<<[]()->bool{return false;}//0
				<<[]()->bool{return false;}//1
				<<[]()->bool{return false;}//2
				<<[]()->bool{return false;}//3
				<<[]()->bool{return false;}//4
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut5",5);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1,0.279,5.446);
						cut->SetPoint(2,0.279,5.217);
						cut->SetPoint(3,0.323,5.217);
						cut->SetPoint(4,0.323,5.446);
						cut->SetPoint(5,0.279,5.446);
					}
					return cut->IsInside(Ek_he(T,P),Th_he(T,P));
				}//5
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut6",5);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1,0.264,5.217);
						cut->SetPoint(2,0.336,5.217);
						cut->SetPoint(3,0.310,5.962);
						cut->SetPoint(4,0.290,5.962);
						cut->SetPoint(5,0.264,5.217);
					}
					return cut->IsInside(Ek_he(T,P),Th_he(T,P));
				}//6
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut7",7);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1,0.263,5.217);
						cut->SetPoint(2,0.342,5.217);
						cut->SetPoint(3,0.332,5.962);
						cut->SetPoint(4,0.311,6.420);
						cut->SetPoint(5,0.286,6.420);
						cut->SetPoint(6,0.273,5.962);
						cut->SetPoint(7,0.263,5.217);
					}
					return cut->IsInside(Ek_he(T,P),Th_he(T,P));
				}//7
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut8",10);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1 ,0.259,5.217);
						cut->SetPoint(2 ,0.282,5.217);
						cut->SetPoint(3 ,0.303,5.618);
						cut->SetPoint(4 ,0.322,5.217);
						cut->SetPoint(5 ,0.354,5.217);
						cut->SetPoint(6 ,0.338,6.134);
						cut->SetPoint(7 ,0.317,6.592);
						cut->SetPoint(8 ,0.291,6.592);
						cut->SetPoint(9 ,0.267,6.134);
						cut->SetPoint(10,0.259,5.217);
					}
					return cut->IsInside(Ek_he(T,P),Th_he(T,P));
				}//8
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut9",10);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1 ,0.252,5.217);
						cut->SetPoint(2 ,0.278,5.217);
						cut->SetPoint(3 ,0.302,5.790);
						cut->SetPoint(4 ,0.329,5.217);
						cut->SetPoint(5 ,0.362,5.217);
						cut->SetPoint(6 ,0.349,6.306);
						cut->SetPoint(7 ,0.320,6.879);
						cut->SetPoint(8 ,0.286,6.879);
						cut->SetPoint(9 ,0.263,6.306);
						cut->SetPoint(10,0.252,5.217);
					}
					return cut->IsInside(Ek_he(T,P),Th_he(T,P));
				}//9
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut10",10);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1 ,0.249,5.217);
						cut->SetPoint(2 ,0.270,5.217);
						cut->SetPoint(3 ,0.304,6.191);
						cut->SetPoint(4 ,0.339,5.217);
						cut->SetPoint(5 ,0.368,5.217);
						cut->SetPoint(6 ,0.353,6.535);
						cut->SetPoint(7 ,0.318,7.166);
						cut->SetPoint(8 ,0.293,7.166);
						cut->SetPoint(9 ,0.263,6.535);
						cut->SetPoint(10,0.249,5.217);
					}
					return cut->IsInside(Ek_he(T,P),Th_he(T,P));
				}//10
				<<[](WTrack&T,const vector<double>&P)->bool{
					static TCutG *cut=nullptr;
					if(cut==nullptr){
						cut=new TCutG("kincut11",11);
						cut->SetVarX("Ek");
						cut->SetVarY("Theta");
						cut->SetPoint(1 ,0.248,5.217);
						cut->SetPoint(2 ,0.268,5.217);
						cut->SetPoint(3 ,0.304,6.306);
						cut->SetPoint(4 ,0.343,5.217);
						cut->SetPoint(5 ,0.379,5.217);
						cut->SetPoint(6 ,0.365,6.535);
						cut->SetPoint(7 ,0.320,7.452);
						cut->SetPoint(8 ,0.292,7.452);
						cut->SetPoint(9 ,0.253,6.535);
						cut->SetPoint(10,0.248,5.217);
					}
					return cut->IsInside(Ek_he(T,P),Th_he(T,P));
				}//11
			)
			<<Forward::Get().CreateMarker(dirname(),"5-Kinematic cut")
			<<make_shared<Hist1D>(dirname(),"5-Kinematic cut",Q);
	}
	shared_ptr<AbstractChain> KinematicHe3Test(const Analysis&data,const Axis&Q,bool MC){
		auto res=make_shared<Chain>()<<make_shared<SetOfHists2D>(dirname(),"Kinematic-reconstructed",Q,Ek_he,Th_he);
		if(MC){
			Axis Ev([&data]()->double{return data.FromFirstVertex(kHe3).E;},Ek_he);
			Axis Tv([&data]()->double{return data.FromFirstVertex(kHe3).Th;},Th_he);
			res<<make_shared<SetOfHists2D>(dirname(),"Kinematic-vertex",Q,Ev,Tv);
		}
		return res;
	}
	///Reaction analysis types visible from reactions.h
	Analysis* He3_forward_analyse(He3Modification mode){
		auto res=Prepare(mode);auto Q=Q_axis(res);
		res->EventProcessing()<<make_shared<Hist1D>(dirname(),"0-Reference",Q);
		res->TrackTypeProcess(kFDC)<<(make_shared<ChainCheck>()
			<<ReconstructionProcess(*res,Q)<<He3Eta_kin_cut(*res,Q)
			<<MissingMass(*res,Q)<<KinematicHe3Test(*res,Q,mode==forEta)
		);
		return res;
	}
	Analysis* He3_forward_reconstruction(He3Modification mode){
		auto res=Prepare(mode);auto Q=Q_axis(res);
		res->EventProcessing()<<make_shared<Hist1D>(dirname(),"0-Reference",Q);
		res->TrackTypeProcess(kFDC)<<(make_shared<ChainCheck>()
			<<ReconstructionProcess(*res,Q)
			<<KinematicHe3Test(*res,Q,mode==forEta)
		);
		return res;
	}
}
