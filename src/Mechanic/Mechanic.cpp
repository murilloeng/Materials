//std
#include <cmath>

//Materials
#include "Materials/inc/Mechanic/Stress.hpp"
#include "Materials/inc/Mechanic/States.hpp"
#include "Materials/inc/Mechanic/Mechanic.hpp"

namespace materials
{
	//constructors
	Mechanic::Mechanic(void)
	{
		return;
	}

	//destructor
	Mechanic::~Mechanic(void)
	{
		return;
	}

	//serialization
	void Mechanic::load(FILE* file)
	{
		Material::load(file);
		if(fscanf(file, "%lf", &m_break_strain) != 1)
		{
			printf("\tError: Unable to load Mechanic material!\n");
			exit(EXIT_FAILURE);
		}
	}
	void Mechanic::save(FILE* file) const
	{
		Material::save(file);
		fprintf(file, "%+.6e ", m_break_strain);
	}

	//properties
	bool Mechanic::fracture(void)
	{
		return m_break;
	}
	bool Mechanic::fracture(bool fracture)
	{
		return m_break = fracture;
	}
	
	bool Mechanic::hardening(void)
	{
		return m_hardening;
	}
	bool Mechanic::hardening(bool hardening)
	{
		return m_hardening = hardening;
	}

	bool Mechanic::inelastic(void)
	{
		return m_inelastic;
	}
	bool Mechanic::inelastic(bool inelastic)
	{
		return m_inelastic = inelastic;
	}

	double Mechanic::break_strain(void) const
	{
		return m_break_strain;
	}
	double Mechanic::break_strain(double break_strain)
	{
		return m_break_strain = break_strain;
	}

	double Mechanic::bulk_modulus(void) const
	{
		const double n = poisson_ratio();
		const double E = elastic_modulus();
		return E / (1 - 2 * n) / 3;
	}
	double Mechanic::shear_modulus(void) const
	{
		const double n = poisson_ratio();
		const double E = elastic_modulus();
		return E / (1 + n) / 2;
	}

	//data
	bool Mechanic::m_break = false;
	bool Mechanic::m_hardening = false;
	bool Mechanic::m_inelastic = false;
}