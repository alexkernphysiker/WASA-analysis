// this file is distributed under 
// MIT license
#include <unistd.h>
#include <math.h>
#include <math_h/interpolate.h>
#include <math_h/error.h>
#include "particles.h"
#include "reactions.h"
#include <experiment_conv.h>
using namespace MathTemplates;
Reaction::Reaction(const Particle& p, const Particle& t, const std::initializer_list< Particle >& products)
:m_projectile(p),m_target(t){
	for(const Particle&item:products)
		m_products.push_back(item);
}
Reaction::Reaction(const Reaction& source)
:m_projectile(source.m_projectile),m_target(source.m_target){
	for(const Particle&item:source.m_products)
		m_products.push_back(item);
}
Reaction::~Reaction(){}
const Particle& Reaction::projectile()const{return m_projectile;}
const Particle& Reaction::target()const{return m_target;}
const std::vector<Particle>& Reaction::products()const{return m_products;}
const double Reaction::M_before() const{return m_projectile.mass()+m_target.mass();}
const double Reaction::M_after() const{
	double res=0;
	for(const Particle&item:m_products)res+=item.mass();
	return res;
}

const double Reaction::EThreshold() const{
	return M_after()-M_before();
}
const double Reaction::PThreshold() const{
	return m_projectile.E2P(EThreshold());
}
const double Reaction::E2Q(const double E) const{
	return sqrt(pow(E+M_before(),2)-pow(m_projectile.E2P(E),2))-M_after();
}
const double Reaction::P2Q(const double P) const{
	return sqrt(pow(m_projectile.P2E(P)+M_before(),2)-pow(P,2))-M_after();
}
const double Reaction::PbEr2Theta(const double Pbeam, const double Ereg) const{
	if(m_products.size()!=2)
		throw Exception<Reaction>("Calculation of theta is available only for binary reactions");
	double E=m_projectile.P2E(Pbeam);
	double P=Pbeam;
	double Q=P2Q(P);
	double beta = P/(E+M_before());
	double gamma = 1./sqrt(1.-pow(beta,2));
	double T_reg_cm = Q/2.*(Q+2*m_products[1].mass())/(Q+M_after());
	double beta_reg_cm = sqrt(pow(T_reg_cm,2)+2*T_reg_cm*m_products[0].mass())/(T_reg_cm+m_products[0].mass());
	double p_reg_cm = beta_reg_cm*(m_products[0].mass()+T_reg_cm);
	double theta_reg_cm = acos((Ereg-(gamma-1)*m_products[0].mass()-gamma*T_reg_cm)/(gamma*beta*p_reg_cm));
	return atan(sin(theta_reg_cm)/(gamma*(cos(theta_reg_cm)+beta/beta_reg_cm)));
}

double sigmaHe3eta(const double Q){
	return 400;
}
double sigmaHe3pi0pi0pi0(const double E){
	return 1000;
}
double sigmaHe3pi0pi0(const double E){
	return 25000;
}
