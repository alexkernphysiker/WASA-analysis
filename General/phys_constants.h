#ifndef AXWBNBYL
#define AXWBNBYL
// Physical constants
const double m_p=0.938272;   //[GeV]
const double m_n=0.93956;    //[GeV]
const double m_d=1.875613;   //[GeV]
const double m_3He=2.808950; //[GeV]
const double m_eta=0.5478; //[GeV]
const double m_pic=0.1396;  //[GeV]
const double m_pi0=0.1350;  //[GeV]
const int c_n=0;
const int c_H=1;
const int c_He=2;
//Conditions of experiment
#define ALLRUNS int runindex=45873;runindex<=46884;runindex++
const double p_beam_low=1.426;
const double p_beam_hi=1.635;
const double p_he3_eta_threshold=1.5727;
const unsigned int beam_momenta_bins=30;
//calsulational designations
inline double NormPhi(double p){
	const double twopi=2*3.1415926;
	double phi=p;
	while(phi<0)phi+=twopi;
	while(phi>=twopi)phi-=twopi;
	return phi;
}
#endif