// this file is distributed under 
// MIT license
#ifndef ______PARTICLES_H_____
#	define ______PARTICLES_H_____
class Particle{
public:
	Particle(const Particle&source);
	virtual ~Particle();
	const double mass()const;
	const int charge()const;
	const double E2P(const double E)const;
	const double P2E(const double P)const;
	static const Particle&n();
	static const Particle&p();
	static const Particle&d();
	static const Particle&he3();
	static const Particle&eta();
	static const Particle&pi0();
	static const Particle&pi_plus();
	static const Particle&pi_minus();
private:
	Particle(double m, int c);
	double m_mass;
	int m_charge;
};

#endif