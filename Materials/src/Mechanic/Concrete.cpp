//std
#include <cmath>

//math
#include "Math/Math/inc/misc/bits.hpp"
#include "Math/Math/inc/linear/vector.hpp"

//materials
#include "Materials/Materials/inc/Types.hpp"
#include "Materials/Materials/inc/Mechanic/Stress.hpp"
#include "Materials/Materials/inc/Mechanic/Concrete.hpp"

namespace materials
{
	//constructors
	Concrete::Concrete(void) :
		m_model(Model::drucker_prager), m_elastic_modulus(2.61e+10), 
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
	Concrete::Model Concrete::model(void) const
	{
		return m_model;
	}
	Concrete::Model Concrete::model(Concrete::Model model)
	{
		return m_model = model;
	}

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
	double Concrete::plastic_stress(double) const
	{
		//data
		const double st = m_yield_stress_tension;
		const double sb = m_yield_stress_biaxial;
		const double sc = m_yield_stress_compression;
		//data
		const double Adp = 2 * st * sc / (st + sc);
		const double Abp = st * sc * sb * (8 * sb - 3 * sc + st) / (st + sc) / (2 * sb + st) / (2 * sb - sc);
		//stress
		return 
			m_model == Model::drucker_prager ? Adp: 
			m_model == Model::bresler_pister ? Abp : 0;
	}
	double Concrete::plastic_gradient(double) const
	{
		return 0;
	}

	//yield criterion
	double Concrete::yield_function(const double* s, uint32_t st) const
	{
		//data
		const double sa = m_yield_stress_tension;
		const double sb = m_yield_stress_biaxial;
		const double sc = m_yield_stress_compression;
		const double I1 = Stress(s, st).invariant_I1();
		const double se = Stress(s, st).von_mises_stress();
		//data
		const double Bdp = -(sc - sa) / (sa + sc);
		const double Cbp = +(3 * sb * sa - 2 * sc * sa - sb * sc) / (sa + sc) / (2 * sb + sa) / (2 * sb - sc);
		const double Bbp = +(sc - sa) * (sb * sc + sb * sa - sc * sa - 4 * sb * sb) / (sa + sc) / (2 * sb + sa) / (2 * sb - sc);
		//yield stress
		return 
			m_model == Model::drucker_prager ? se - Bdp * I1 : 
			m_model == Model::bresler_pister ? se - Bbp * I1 - Cbp * I1 * I1 : 0;
	}
	void Concrete::yield_hessian(double* ddf, const double* s, uint32_t st) const
	{
		//data
		double dI1[6];
		const uint32_t ns = math::bit_count(st);
		const double sa = m_yield_stress_tension;
		const double sb = m_yield_stress_biaxial;
		const double sc = m_yield_stress_compression;
		const double Cbp = +(3 * sb * sa - 2 * sc * sa - sb * sc) / (sa + sc) / (2 * sb + sa) / (2 * sb - sc);
		//hessian
		Stress(s, st).von_mises_hessian(ddf);
		math::matrix(ddf, ns, ns) -= Cbp * math::vector(dI1, ns).outer(math::vector(dI1, ns));
	}
	void Concrete::yield_gradient(double* df, const double* s, uint32_t st) const
	{
		//data
		double dI1[6];
		const uint32_t ns = math::bit_count(st);
		const double sa = m_yield_stress_tension;
		const double sb = m_yield_stress_biaxial;
		const double sc = m_yield_stress_compression;
		const double I1 = Stress(s, st).invariant_I1();
		//data
		const double Bdp = (sa - sc) / (sa + sc);
		const double Cbp = +(3 * sb * sa - 2 * sc * sa - sb * sc) / (sa + sc) / (2 * sb + sa) / (2 * sb - sc);
		const double Bbp = +(sc - sa) * (sb * sc + sb * sa - sc * sa - 4 * sb * sb) / (sa + sc) / (2 * sb + sa) / (2 * sb - sc);
		//yield gradient
		Stress(s, st).gradient_I1(dI1);
		Stress(s, st).von_mises_gradient(df);
		if(m_model == Model::drucker_prager)
		{
			math::vector(df, ns) -= Bdp * math::vector(dI1, ns);
		}
		if(m_model == Model::bresler_pister)
		{
			math::vector(df, ns) -= (Bbp + 2 * I1 * Cbp) * math::vector(dI1, ns);
		}
	}
}