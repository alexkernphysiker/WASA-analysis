// this file is distributed under 
// MIT license
#ifndef AXWBNBYL
#define AXWBNBYL

#define ALLRUNS int runindex=45873;runindex<=46884;runindex++
const double p_beam_low=1.426;
const double p_beam_hi=1.635;
//calsulational designations
inline double NormPhi(double p){
	const double twopi=2*3.1415926;
	double phi=p;
	while(phi<0)phi+=twopi;
	while(phi>=twopi)phi-=twopi;
	return phi;
}
#endif