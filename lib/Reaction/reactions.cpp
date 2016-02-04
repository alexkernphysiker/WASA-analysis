// this file is distributed under 
// MIT license
#include <unistd.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "particles.h"
#include "reactions.h"
double Q_He3eta(double pBeam){
	TVector3 p_beam;
	p_beam.SetMagThetaPhi(pBeam,0,0);
	TLorentzVector P_Beam;
	TLorentzVector P_Target;
	P_Beam.SetVectM(p_beam,Particle::p().mass_GeV());
	TVector3 ptarget;
	ptarget.SetMagThetaPhi(0,0,0);
	P_Target.SetVectM(ptarget,Particle::d().mass_GeV());
	TLorentzVector P_Total=P_Beam+P_Target;
	return P_Total.M()-(Particle::he3().mass_GeV()+Particle::eta().mass_GeV());
}
double sigmaHe3eta(double q){
	if(q<=0)return 0;
	else return 400;
}