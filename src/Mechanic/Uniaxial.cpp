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

	//return mapping
	void Uniaxial::stress(double& s, double e, const Point& mp) const
	{
		//data
		const double E = m_elastic_modulus;
		const double ep = mp.m_plastic_strain_new[0];
		//stress
		s = !m_inelastic ? E * e : 
			mp.m_state_new == state::broken ? 0 :
			mp.m_state_new == state::elastic ? E * e :
			mp.m_state_new == state::plastic ? E * (e - ep) : 0;
	}
	void Uniaxial::stiffness(double& C, double e, const Point& mp) const
	{
		//data
		const double E = m_elastic_modulus;
		const double K = m_plastic_modulus;
		//stiffness
		C = !m_inelastic ? E :
			mp.m_state_new == state::broken ? 0 :
			mp.m_state_new == state::elastic ? E :
			mp.m_state_new == state::plastic ? E * K / (E + K) : 0;
	}
	void Uniaxial::return_mapping(double& s, double& C, double e, Point& mp) const
	{
		//data
		const double sy = m_yield_stress;
		const double sb = m_buckling_stress;
		const double E0 = m_elastic_modulus;
		const double Kp = m_plastic_modulus;
		//check
		if(!m_inelastic && !m_buckling)
		{
			C = E0;
			s = E0 * e;
			return;
		}
		//buckling
		if(m_buckling && fabs(s) > sb)
		{
			C = 0;
			s = math::sign(s) * sb;
			return;
		}
		//inelastic
		if(m_inelastic)
		{
			
		}
	}
}