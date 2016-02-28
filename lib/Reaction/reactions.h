// this file is distributed under 
// MIT license
#ifndef xylgjnjy
#	define xylgjnjy
#include <list>
#include <vector>
#include "particles.h"
const double p_he3_eta_threshold=1.5727;
class Reaction{
public:
	Reaction(const Particle&p,const Particle&t,const std::initializer_list<Particle>&products);
	Reaction(const Reaction&source);
	virtual ~Reaction();
	const Particle&target()const;
	const Particle&projectile()const;
	const std::vector<Particle>&products()const;
	const double M_before()const;
	const double M_after()const;
	const double PThreshold()const;
	const double EThreshold()const;
	const double E2Q(const double E)const;
	const double P2Q(const double P)const;
	const double PbEr2Theta(const double Pbeam,const double Ereg)const;
private:
	Particle m_projectile,m_target;
	std::vector<Particle> m_products;
};
//known cross sections
double sigmaHe3eta(const double Q);
double sigmaHe3pi0pi0pi0(const double E);
double sigmaHe3pi0pi0(const double E);
#endif 