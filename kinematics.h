// this file is distributed under 
// GPL v 3.0 license
#ifndef rtyaxbvk
#define rtyaxbvk
#include <TVector3.h>
#include <TLorentzVector.h>
#include "phys_constants.h"
inline double Q_He3eta(double pBeam){
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
#endif 