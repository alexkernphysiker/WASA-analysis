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
	shared_ptr<ICalculation> FRH1_cut(const Analysis&data){
		auto res=make_shared<CalcChain>();
		res<<[](WTrack&T,Results&)->bool{
			return Forward::Get().StoppingLayer(T)==kFRH1;
		};
		res<<[](WTrack&T,Results&)->bool{
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
		};
		res<<Forward::Get().CreateMarker(dirname(),"2.5-FRH1");
		res<<[&data](WTrack&track,Results&P)->bool{
			//Achtung - static
			static FitBasedReconstruction<Reconstruction::He3EnergyFRH1,WTrack&> energy(
				"He3.E.FRH1",
				{
					[](WTrack&track){return Forward::Get()[kFRH1].Edep(track);},
					[](WTrack&track){return track.Theta();}
				},
				[&data](WTrack&){return data.FromFirstVertex(kHe3).E;}
			);
			P.push_back(energy.Reconstruct(track));
			return true;
		};
		return res;
	}
	shared_ptr<ICalculation> FRH2_cut(const Analysis&data){
		auto res=make_shared<CalcChain>();
		res<<[](WTrack&T,Results&)->bool{
			return Forward::Get().StoppingLayer(T)==kFRH2;
		};
		res<<[](WTrack&T,Results&)->bool{
			vector<double> Ed={
				Forward::Get()[kFRH1].Edep(T),
				Forward::Get()[kFRH2].Edep(T)
			};
			return (Ed[0]>(0.25-0.417*Ed[1]))&&(Ed[0]<(0.35-0.417*Ed[1]))&&(Ed[1]<0.22);
		};
		res<<Forward::Get().CreateMarker(dirname(),"2.5-FRH2");
		res<<[&data](WTrack&track,Results&P)->bool{
			//Achtung - static
			static FitBasedReconstruction<Reconstruction::He3EnergyFRH2,WTrack&> energy(
				"He3.E.FRH2",
				{
					[](WTrack&T){return Forward::Get()[kFRH1].Edep(T)+Forward::Get()[kFRH2].Edep(T);},
					[](WTrack&T){return T.Theta();}
				},
				[&data](WTrack&){return data.FromFirstVertex(kHe3).E;}
			);
			P.push_back(energy.Reconstruct(track));
			return true;
		};
		return res;
	}
	shared_ptr<ICalculation> ReconstructionProcess(const Analysis&data,const Axis<>&Q){
		auto res=make_shared<CalcChain>();
		res<<Forward::Get().CreateMarker(dirname(),"1-AllTracks");
		res<<make_shared<Hist1D<>>(dirname(),"1-AllTracks",Q);
			
		res<<[](WTrack&T,Results&)->bool{
			//ToDo: replace by more reasonable condition
			return (T.Theta()!=0.125);
		};
			
		res<<Forward::Get().CreateMarker(dirname(),"2-FPC");
		res<<make_shared<Hist1D<>>(dirname(),"2-FPC",Q);
			
		auto cuts=make_shared<CalcOr>();
		cuts<<FRH1_cut(data)<<FRH2_cut(data);
		res<<cuts;
			
		res<<Forward::Get().CreateMarker(dirname(),"3-AllCuts");
		res<<make_shared<Hist1D<>>(dirname(),"3-AllCuts",Q);
			
		res<<[](WTrack&T,Results&R)->bool{
			R.push_back(T.Theta());
			return true;
		};
		res<<[](WTrack&T,Results&R)->bool{
			R.push_back(T.Phi());
			return true;
		};
		res<<[](WTrack&,Results&R)->bool{
			return isfinite(R[0])&&isfinite(R[1])&&isfinite(R[2]);
		};
			
		res<<Forward::Get().CreateMarker(dirname(),"4-Reconstructed");
		res<<make_shared<Hist1D<>>(dirname(),"4-Reconstructed",Q);
		return res;
	}
	shared_ptr<ICalculation> MissingMass(const Analysis&data,const Axis<>&Q){
		Axis<WTrack&,Results&> Q_([Q](WTrack&,Results&){return Q();},Q);
		Axis<WTrack&,Results&> M([](WTrack&,Results&P)->double{return P[3];},m_pi0-0.5,m_eta+0.5,100);
		auto res=make_shared<CalcChain>();
		res<<[&data](WTrack&,Results&P)->bool{
			double Ekin=P[0];
			double theta=P[1];
			double phi=P[2];
			double p=sqrt(Ekin*(Ekin+2*m_3He));
			TVector3 p_He3;
			p_He3.SetMagThetaPhi(p,theta,phi);
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
			P.push_back(P_Missing.M());
			return true;
		};
		res<<make_shared<SetOfHists1D<WTrack&,Results&>>(dirname(),"MissingMass",Q_,M);
		return res;
	}
	shared_ptr<ICalculation> He3Eta_cut(){
		auto res=make_shared<CalcChain>();
		res<<[](WTrack&track,Results&)->bool{
			if(Forward::Get().StoppingLayer(track)!=kFRH1)return false;
			double E=Forward::Get()[kFRH1].Edep(track);
			return (E>0.08)&&(E<0.22);
		};
		return res;
	}
	shared_ptr<IResAnalysis> KinematicHe3Test(const Analysis&data,He3Modification mode){
		Axis<const Results&> Bkin([&data](const Results&)->double{return 1000.0*Q_He3eta(data.PBeam());},Bkin);
		Axis<const Results&> Ev([&data](const Results&)->double{return data.FromFirstVertex(kHe3).E;},0.1,0.6,50);
		Axis<const Results&> Tv([&data](const Results&)->double{return data.FromFirstVertex(kHe3).Th;},0.10,0.16,50);
		Axis<const Results&> Ek([](const Results&P)->double{return P[0];},Ev);
		Axis<const Results&> Th([](const Results&P)->double{return P[1];},Tv);
		auto res=make_shared<AnalysisSet>();
		res<<make_shared<SetOfHists2D<const Results&>>(dirname(),"Kinematic-reconstructed",Bkin,Ek,Th);
		if(mode==forEta)
			res<<make_shared<SetOfHists2D<const Results&>>(dirname(),"Kinematic-vertex",Bkin,Ev,Tv);
		return res;
	}
	
	///Reaction analysis types visible from reactions.h
	Analysis* He3_forward_analyse(He3Modification mode){
		auto res=Prepare(mode);
		Axis<> Q([res]()->double{return 1000.0*Q_He3eta(res->PBeam());},0.0,30.0,12);
		auto chain=make_shared<TrackCalcChain>();
		chain<<ReconstructionProcess(*res,Q)<<He3Eta_cut();
		chain<<KinematicHe3Test(*res,mode)<<MissingMass(*res,Q);
		res->TrackTypeProcess(kFDC).Add(chain);
		return res;
	}
	Analysis* He3_forward_reconstruction(He3Modification mode){
		auto res=Prepare(mode);
		Axis<> Q([res]()->double{return 1000.0*Q_He3eta(res->PBeam());},0.0,30.0,12);
		auto chain=make_shared<TrackCalcChain>();
		chain<<ReconstructionProcess(*res,Q)<<KinematicHe3Test(*res,mode);
		res->TrackTypeProcess(kFDC).Add(chain);
		return res;
	}
}
