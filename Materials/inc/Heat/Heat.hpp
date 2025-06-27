#pragma once

//materials
#include "Materials/Materials/inc/Material.hpp"

namespace materials
{
	class Heat : public Material
	{
		friend class Material;

	protected:
		//constructors
		Heat(void);

		//destructor
		~Heat(void) override;

		//serialization
		void load(FILE*) override;
		void save(FILE*) const override;

	public:
		//data
		double capacity(double);
		double capacity(void) const;
		
		double condutivity(double);
		double condutivity(void) const;

		//type
		materials::type type(void) const override;

	protected:
		//data
		double m_capacity;
		double m_condutivity;
	};
}