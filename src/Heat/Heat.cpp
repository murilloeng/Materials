//std
#include <stdexcept>

//Materials
#include "Materials/inc/Heat/Heat.hpp"

namespace materials
{
	//constructors
	Heat::Heat(void) : m_capacity(1), m_condutivity(1)
	{
		return;
	}

	//destructor
	Heat::~Heat(void)
	{
		return;
	}

	//serialization
	void Heat::load(FILE* file)
	{
		Material::load(file);
		if(fscanf(file, "%lf %lf", &m_condutivity, &m_capacity) != 2)
		{
			throw std::runtime_error("Error: Unable to load Heat material!");
		}
	}
	void Heat::save(FILE* file) const
	{
		Material::save(file);
		if(fprintf(file, "%+.6e %+.6e ", m_condutivity, m_capacity) < 0)
		{
			throw std::runtime_error("Error: Unable to save Heat material!");
		}
	}

	//data
	double Heat::capacity(void) const
	{
		return m_capacity;
	}
	double Heat::capacity(double capacity)
	{
		return m_capacity = capacity;
	}
	double Heat::condutivity(void) const
	{
		return m_condutivity;
	}
	double Heat::condutivity(double condutivity)
	{
		return m_condutivity = condutivity;
	}
}