//std
#include <stdexcept>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Miscellaneous/bits.hpp"

//Materials
#include "Materials/inc/Mechanic/Steel.hpp"
#include "Materials/inc/Mechanic/Stress.hpp"

namespace materials
{
	//constructors
	Steel::Steel(void) :
		m_yield_stress(400e6), m_poisson_ratio(0.3), m_plastic_modulus(0), m_elastic_modulus(200e9)
	{
		m_break_strain = 2.00e-01;
		m_specific_mass = 7.85e+03;
	}

	//destructor
	Steel::~Steel(void)
	{
		return;
	}

	//serialization
	void Steel::load(FILE* file)
	{
		Mechanic::load(file);
		const char* format = "%lf %lf %lf %lf";
		if(fscanf(file, format, &m_elastic_modulus, &m_poisson_ratio, &m_yield_stress, &m_plastic_modulus) != 4)
		{
			throw std::runtime_error("Error: Unable to load Steel material!");
		}
	}
	void Steel::save(FILE* file) const
	{
		Mechanic::save(file);
		const char* format = "%+.6e %+.6e %+.6e %+.6e ";
		if(fprintf(file, format, m_elastic_modulus, m_poisson_ratio, m_yield_stress, m_plastic_modulus) < 0)
		{
			throw std::runtime_error("Error: Unable to save Steel material!");
		}
	}

	//data
	double Steel::yield_stress(void) const
	{
		return m_yield_stress;
	}
	double Steel::yield_stress(double yield_stress)
	{
		return m_yield_stress = yield_stress;
	}

	double Steel::plastic_modulus(void) const
	{
		return m_plastic_modulus;
	}
	double Steel::plastic_modulus(double plastic_modulus)
	{
		return m_plastic_modulus = plastic_modulus;
	}

	double Steel::poisson_ratio(void) const
	{
		return m_poisson_ratio;
	}
	double Steel::poisson_ratio(double poisson_ratio)
	{
		return m_poisson_ratio = poisson_ratio;
	}

	double Steel::elastic_modulus(void) const
	{
		return m_elastic_modulus;
	}
	double Steel::elastic_modulus(double elastic_modulus)
	{
		return m_elastic_modulus = elastic_modulus;
	}

	double Steel::thermal_coefficient(void) const
	{
		return 1.2e-5;
	}

	//plastic stress
	double Steel::plastic_stress(double p) const
	{
		return m_yield_stress + m_plastic_modulus * p;
	}
	double Steel::plastic_gradient(double p) const
	{
		return m_plastic_modulus;
	}

	//yield criterion
	double Steel::yield_function(const double* s, uint32_t st) const
	{
		return Stress(s, st).von_mises_stress();
	}
	void Steel::yield_gradient(double* df, const double* s, uint32_t st) const
	{
		Stress(s, st).von_mises_gradient(df);
	}
	void Steel::yield_hessian(double* ddf, const double* s, uint32_t st) const
	{
		Stress(s, st).von_mises_hessian(ddf);
	}
}