//std
#include <cmath>

//Math
#include "Math/inc/Miscellaneous/util.hpp"

//Materials
#include "Materials/inc/Mechanic/States.hpp"
#include "Materials/inc/Mechanic/Uniaxial.hpp"

namespace materials
{
	//constructors
	Uniaxial::Uniaxial(void)
	{
		return;
	}

	//destructor
	Uniaxial::~Uniaxial(void)
	{
		return;
	}

	//data
	bool Uniaxial::buckling(void)
	{
		return m_buckling;
	}
	bool Uniaxial::buckling(bool buckling)
	{
		return m_buckling = buckling;
	}

	double Uniaxial::yield_stress(void) const
	{
		return m_yield_stress;
	}
	double Uniaxial::yield_stress(double yield_stress)
	{
		return m_yield_stress;
	}

	double Uniaxial::buckling_stress(void) const
	{
		return m_buckling_stress;
	}
	double Uniaxial::buckling_stress(double buckling_stress)
	{
		return m_buckling_stress = buckling_stress;
	}

	double Uniaxial::elastic_modulus(void) const
	{
		return m_elastic_modulus;
	}
	double Uniaxial::elastic_modulus(double elastic_modulus)
	{
		return m_elastic_modulus = elastic_modulus;
	}

	double Uniaxial::plastic_modulus(void) const
	{
		return m_plastic_modulus;
	}
	double Uniaxial::plastic_modulus(double plastic_modulus)
	{
		return m_plastic_modulus = plastic_modulus;
	}

	Uniaxial::Function Uniaxial::yield_stress_function(void) const
	{
		return m_yield_stress_function;
	}
	Uniaxial::Function Uniaxial::yield_stress_function(Function yield_stress_function)
	{
		return m_yield_stress_function = yield_stress_function;
	}

	Uniaxial::Function Uniaxial::plastic_modulus_function(void) const
	{
		return m_plastic_modulus_function;
	}
	Uniaxial::Function Uniaxial::plastic_modulus_function(Function plastic_modulus_function)
	{
		return m_plastic_modulus_function = plastic_modulus_function;
	}

	//return mapping
	void Uniaxial::stress(double& s, double e, const Point& point) const
	{
		//data
		const double E = m_elastic_modulus;
		const double ep = point.m_plastic_strain_new[0];
		//stress
		s = !m_inelastic ? E * e : 
			point.m_state_new == state::broken ? 0 : E * (e - ep);
	}
	void Uniaxial::stiffness(double& C, double e, const Point& point) const
	{
		//data
		const double E = m_elastic_modulus;
		const double p = point.m_hardening_new;
		const double K = compute_plastic_modulus(p);
		//stiffness
		C = !m_inelastic ? E :
			point.m_state_new == state::broken ? 0 :
			point.m_state_new == state::elastic ? E : E * K / (E + K);
	}
	void Uniaxial::return_mapping(double& s, double& C, double e, Point& mp) const
	{
		//data
		double f, sy, Kp;
		const double E = m_elastic_modulus;
		const double sr = compute_yield_stress(0);
		//new state
		double& p_new = mp.m_hardening_new;
		double& ep_new = mp.m_plastic_strain_new[0];
		//old state
		const double p_old = mp.m_hardening_old;
		const double ep_old = mp.m_plastic_strain_old[0];
		//predictor
		C = E;
		p_new = p_old;
		ep_new = ep_old;
		s = E * (e - ep_new);
		if(fabs(s) < compute_yield_stress(p_new)) return;
		//corrector
		for(uint32_t i = 0; i < 10; i++)
		{
			//compute
			sy = compute_yield_stress(p_new);
			Kp = compute_plastic_modulus(p_new);
			ep_new = ep_old + (p_new - p_old) * math::sign(s);
			//check
			f = fabs(E * (e - ep_new)) - sy;
			if(fabs(f) < 1e-5 * sr) break;
			//update
			p_new += f / (E + Kp);
		}
	}

	//compute
	double Uniaxial::compute_yield_stress(double p) const
	{
		return m_yield_stress_function ? m_yield_stress_function(p) : m_yield_stress;
	}
	double Uniaxial::compute_plastic_modulus(double p) const
	{
		return m_plastic_modulus_function ? m_plastic_modulus_function(p) : m_plastic_modulus;
	}

	bool Uniaxial::m_buckling = true;
}