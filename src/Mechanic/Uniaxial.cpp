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
	Uniaxial::Uniaxial(void) : 
		m_yield_stress{0}, m_buckling_stress{0}, m_elastic_modulus{0}, m_plastic_modulus{0}, m_thermal_coefficient{0}
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
		return m_yield_stress = yield_stress;
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

	double Uniaxial::thermal_coefficient(void) const
	{
		return m_thermal_coefficient;
	}
	double Uniaxial::thermal_coefficient(double thermal_coefficient)
	{
		return m_thermal_coefficient = thermal_coefficient;
	}

	double Uniaxial::plastic_stress(double p) const
	{
		return m_yield_stress_function ? m_yield_stress_function(p) : m_yield_stress + m_plastic_modulus * p;
	}
	double Uniaxial::plastic_gradient(double p) const
	{
		return m_plastic_modulus_function ? m_plastic_modulus_function(p) : m_plastic_modulus;
	}

	double Uniaxial::poisson_ratio(void) const
	{
		return 0;
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
		const double K = plastic_gradient(point.m_hardening_new);
		//stiffness
		C = !m_inelastic ? E :
			point.m_state_new == state::broken ? 0 :
			point.m_state_new == state::elastic ? E : E * K / (E + K);
	}
	bool Uniaxial::return_mapping(double& s, double& C, double e, double sr, Point& point) const
	{
		return 
			return_mapping_elastic(s, C, e, sr, point) ? true :
			return_mapping_fracture(s, C, e, sr, point) ? true :
			return_mapping_predictor(s, C, e, sr, point) ? true : 
			return_mapping_corrector(s, C, e, sr, point);
	}

	//return mapping
	bool Uniaxial::return_mapping_elastic(double& s, double& C, double e, double sr, Point& point) const
	{
		//data
		C = m_elastic_modulus;
		s = m_elastic_modulus * e + sr;
		point.m_state_new = state::elastic;
		//return
		return !m_inelastic;
	}
	bool Uniaxial::return_mapping_buckling(double& s, double& C, double e, double sr, Point& point) const
	{
		return false;
	}
	bool Uniaxial::return_mapping_fracture(double& s, double& C, double e, double sr, Point& point) const
	{
		//data
		s = C = 0;
		point.m_state_new = state::broken;
		//return
		return m_break && fabs(e) > m_break_strain;
	}
	bool Uniaxial::return_mapping_predictor(double& s, double& C, double e, double sr, Point& point) const
	{
		//data
		const double E = C = m_elastic_modulus;
		//new state
		double& p_new = point.m_hardening_new;
		double& ep_new = point.m_plastic_strain_new[0];
		//old state
		const double p_old = point.m_hardening_old;
		const double ep_old = point.m_plastic_strain_old[0];
		//predictor
		p_new = p_old;
		ep_new = ep_old;
		s = E * (e - ep_old) + sr;
		point.m_state_new = state::elastic;
		//return
		return fabs(s) < plastic_stress(p_new);
	}
	bool Uniaxial::return_mapping_corrector(double& s, double& C, double e, double sr, Point& point) const
	{
		//data
		double f, sy, Kp;
		const double E = m_elastic_modulus;
		const uint32_t iteration_max = 100;
		const double sa = plastic_stress(0);
		//new state
		double& p_new = point.m_hardening_new;
		double& ep_new = point.m_plastic_strain_new[0];
		//old state
		const double p_old = point.m_hardening_old;
		const double ep_old = point.m_plastic_strain_old[0];
		//corrector
		point.m_state_new = state::plastic;
		for(uint32_t i = 0; i < iteration_max; i++)
		{
			//compute
			sy = plastic_stress(p_new);
			Kp = plastic_gradient(p_new);
			ep_new = ep_old + (p_new - p_old) * math::sign(s);
			//state
			C = E * Kp / (E + Kp);
			s = E * (e - ep_new) + sr;
			//check
			f = fabs(s) - sy;
			if(fabs(f) < 1e-5 * sa) return true;
			//update
			p_new += f / (E + Kp);
		}
		//return
		return false;
	}

	bool Uniaxial::m_buckling = false;
}