// this file is distributed under 
// MIT license
#include <unistd.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "../../phys_constants.h"
#include "he3.h"
namespace Theory{
	double Q_He3eta(double pBeam){
		TVector3 p_beam;
		p_beam.SetMagThetaPhi(pBeam,0,0);
		TLorentzVector P_Beam;
		TLorentzVector P_Target;
		P_Beam.SetVectM(p_beam,m_p);
		TVector3 ptarget;
		ptarget.SetMagThetaPhi(0,0,0);
		P_Target.SetVectM(ptarget,m_d);
		TLorentzVector P_Total=P_Beam+P_Target;
		return P_Total.M()-(m_3He+m_eta);
	}
	double sigmaHe3eta(double p_beam){
		if(p_beam<=p_he3_eta_threshold)return 0;
		else return 400;
	}
};