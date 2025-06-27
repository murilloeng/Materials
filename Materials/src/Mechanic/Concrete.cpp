//std
#include <cmath>

//math
#include "Math/Math/inc/misc/bits.hpp"
#include "Math/Math/inc/linear/vector.hpp"

//materials
#include "Materials/Materials/inc/Types.hpp"
#include "Materials/Materials/inc/Mechanic/stress.hpp"
#include "Materials/Materials/inc/Mechanic/Concrete.hpp"

namespace materials
{
	//constructors
	Concrete::Concrete(void) :
		m_elastic_modulus(2.61e+10), 
		m_yield_stress_tension(3.00e+06), m_yield_stress_compression(3.00e+07)
	{
		m_break_strain = 3.00e-03;
		m_specific_mass = 2.50e+03;
	}

	//destructor
	Concrete::~Concrete(void)
	{
		return;
	}

	//serialization
	void Concrete::load(FILE* file)
	{
		Mechanic::load(file);
		if(fscanf(file, "%lf %lf %lf",
			&m_elastic_modulus, &m_yield_stress_tension, &m_yield_stress_compression) != 3)
		{
			printf("\tError: Unable to load Concrete material!\n");
			exit(EXIT_FAILURE);
		}
	}
	void Concrete::save(FILE* file) const
	{
		Mechanic::save(file);
		fprintf(file, "%+.6e %+.6e %+.6e ",
			m_elastic_modulus,
			m_yield_stress_tension,
			m_yield_stress_compression
		);
	}

	//type
	materials::type Concrete::type(void) const
	{
		return materials::type::concrete;
	}

	//data
	double Concrete::poisson_ratio(void) const
	{
		return 0.2;
	}

	double Concrete::elastic_modulus(double E)
	{
		return m_elastic_modulus = E;
	}
	double Concrete::elastic_modulus(void) const
	{
		return m_elastic_modulus;
	}

	double Concrete::yield_stress_tension(double st)
	{
		return m_yield_stress_tension = st;
	}
	double Concrete::yield_stress_tension(void) const
	{
		return m_yield_stress_tension;
	}

	double Concrete::yield_stress_compression(double sc)
	{
		return m_yield_stress_compression = sc;
	}
	double Concrete::yield_stress_compression(void) const
	{
		return m_yield_stress_compression;
	}

	double Concrete::thermal_coefficient(void) const
	{
		return 1e-5;
	}

	//plastic stress
	double Concrete::plastic_stress(double p) const
	{
		//data
		const double sa = m_yield_stress_tension;
		const double sc = m_yield_stress_compression;
		//stress
		return 2 * sa * sc / (sa + sc);
	}
	double Concrete::plastic_gradient(double p) const
	{
		return 0;
	}

	//yield criterion
	double Concrete::yield_function(const double* s, uint32_t st) const
	{
		//data
		const double sa = m_yield_stress_tension;
		const double sc = m_yield_stress_compression;
		//stress
		const double se = von_mises_stress(s, st);
		const double I1 = stress_invariant_I1(s, st);
		const double B = (sa - sc) / (sa + sc);
		//yield stress
		return se - B * I1;
	}
	void Concrete::yield_hessian(double* dfss, const double* s, uint32_t st) const
	{
		von_mises_hessian(dfss, s, st);
	}
	void Concrete::yield_gradient(double* dfs, const double* s, uint32_t st) const
	{
		//data
		double dI1[6];
		const uint32_t ns = math::bit_count(st);
		const double sa = m_yield_stress_tension;
		const double sc = m_yield_stress_compression;
		//yield gradient
		stress_gradient_I1(dI1, s, st);
		von_mises_gradient(dfs, s, st);
		const double B = (sa - sc) / (sa + sc);
		math::vector(dfs, ns) -= B * math::vector(dI1, ns);
	}
}