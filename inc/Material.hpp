#pragma once

//std
#include <cstdio>
#include <cstdint>

namespace materials
{
	enum class Type : uint32_t;
}

namespace materials
{
	class Material
	{
	protected:
		//constructors
		Material(void);

		//destructor
		virtual ~Material(void);

		//serialization
		virtual void load(FILE*);
		virtual void save(FILE*) const;

	public:
		//data
		double specific_mass(double);
		double specific_mass(void) const;

		//analysis
		virtual void check(void) const;

	protected:
		//data
		double m_specific_mass;
	};
}