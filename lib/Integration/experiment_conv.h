// this file is distributed under 
// MIT license
#ifndef AXWBNBYL
#	define AXWBNBYL
#define ALLRUNS int runindex=45873;runindex<=46884;runindex++
const double p_beam_low=1.426;
const double p_beam_hi=1.635;
const double pbeam_measurement_offset=0.006;
struct trigger{unsigned char number; unsigned long scaling;};
const trigger trigger_he3_forward={.number=10,.scaling=1};

//calsulational designations
inline double NormPhi(double p){
	const double twopi=2*3.1415926;
	double phi=p;
	while(phi<0)phi+=twopi;
	while(phi>=twopi)phi-=twopi;
	return phi;
}
#endif