#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <phys_constants.h>
using namespace std;
#define ParR(x) static_cast<ParamSet&&>(x)
int main(int,char**){
	#include "env.cc"
	size_t n=2;
	double P[]={1.426,1.635};
	for(size_t i=0;i<n;i++){
		TVector3 p_beam;
		p_beam.SetMagThetaPhi(P[i],0,0);
		TLorentzVector P_Beam;
		TLorentzVector P_Target;
		P_Beam.SetVectM(p_beam,m_p);
		TVector3 ptarget;
		ptarget.SetMagThetaPhi(0,0,0);
		P_Target.SetVectM(ptarget,m_d);
		TLorentzVector P_Total=P_Beam+P_Target;
		double Q=P_Total.M()-(m_3He+m_eta);
		printf("P=%f GeV/c => Q=%f GeV\n",P[i],Q);
	}
}