//std
#include <stdexcept>

//Materials
#include "Materials/inc/Material.hpp"

namespace materials
{
	//constructors
	Material::Material(void) : m_specific_mass{0}
	{
		return;
	}

	//destructor
	Material::~Material(void)
	{
		return;
	}

	//serialization
	void Material::load(FILE* file)
	{
		if(fscanf(file, "%lf", &m_specific_mass) != 1)
		{
			throw std::runtime_error("Error: Unable to load Material!");
		}
	}
	void Material::save(FILE* file) const
	{
		if(fprintf(file, "%+.6e ", m_specific_mass) < 0)
		{
			throw std::runtime_error("Error: Unable to save Material!");
		}
	}

	//data
	double Material::specific_mass(void) const
	{
		return m_specific_mass;
	}
	double Material::specific_mass(double specific_mass)
	{
		return m_specific_mass = specific_mass;
	}

	//analysis
	void Material::check(void) const
	{
		if(m_specific_mass < 0)
		{
			throw std::runtime_error("Error: Material with negative specific mass!");
		}
	}
}