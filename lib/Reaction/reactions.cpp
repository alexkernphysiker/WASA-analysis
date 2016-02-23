// this file is distributed under 
// MIT license
#include <unistd.h>
#include <math.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <math_h/interpolate.h>
#include "particles.h"
#include "reactions.h"
#include <experiment_conv.h>
using namespace MathTemplates;
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
double PBeam_He3eta(double Q){
	static LinearInterpolation<double> P_Q=LinearInterpolation<double>(
		Q_He3eta,ChainWithStep(0.0,0.001,3.0)
	).Transponate();
	return P_Q(Q);
}
double Ekin2Theta_He3eta(double Ekin,double p_beam){
	double mbeam=Particle::p().mass_GeV();
	double mtarget=Particle::d().mass_GeV();
	double Eb = sqrt(pow(p_beam,2)+pow(mbeam,2))-mbeam;
	double beta = p_beam/(Eb+mbeam+mtarget);
	double gamma = 1./sqrt(1.-pow(beta,2));
	double Q = Q_He3eta(p_beam);
	double meta=Particle::eta().mass_GeV();
	double mHe3=Particle::he3().mass_GeV();
	double THe3_cm = Q/2.*(Q+2*meta)/(Q+meta+mHe3);
	double betaHe3_cm = sqrt(pow(THe3_cm,2)+2*THe3_cm*mHe3)/(THe3_cm+mHe3);
	double pHe3_cm = betaHe3_cm*(mHe3+THe3_cm);
	double theta_cm = acos((Ekin-(gamma-1)*mHe3-gamma*THe3_cm)/(gamma*beta*pHe3_cm));
	return atan(sin(theta_cm)/(gamma*(cos(theta_cm)+beta/betaHe3_cm)))*180./3.1415926;
}
double sigmaHe3eta(double p_beam){
	//From proposal
	if(p_beam<=p_he3_eta_threshold)return 0;
	else return 400;
}
double sigmaHe3pi0pi0pi0(double p_beam){
	//From proposal
	return 27;
}
double sigmaHe3pi0pi0(double p_beam){
	//From proposal
	return 2800;
}
