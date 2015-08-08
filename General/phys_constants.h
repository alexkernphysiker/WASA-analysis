#ifndef AXWBNBYL
#define AXWBNBYL
const double m_p=0.938272;   //[GeV]
const double m_n=0.93956;    //[GeV]
const double m_d=1.875613;   //[GeV]
const double m_3He=2.808950; //[GeV]
const double m_eta=0.547853; //[GeV]
const double m_pi=0.135;  //[GeV]
const int c_n=0;
const int c_H=1;
const int c_He=2;
const double p_beam_low=1.426;
const double p_beam_hi=1.635;
const double p_he3_eta_threshold=1.5727;
inline double NormPhi(double p){
	double phi=p;
	while(phi<0)phi+=2*3.1415926;
	while(phi>=2*3.1415926)phi-=2*3.1415926;
	return phi;
}

#endif