//std
#include <cstdlib>

//materials
#include "Materials/Materials/inc/Types.hpp"
#include "Materials/Materials/inc/Heat/Heat.hpp"

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
			printf("\tError: Unable to load Heat material!\n");
			exit(EXIT_FAILURE);
		}
	}
	void Heat::save(FILE* file) const
	{
		Material::save(file);
		fprintf(file, "%+.6e %+.6e ", m_condutivity, m_capacity);
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

	//type
	materials::type Heat::type(void) const
	{
		return materials::type::heat;
	}
}