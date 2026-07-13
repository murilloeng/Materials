#pragma once

//Materials
#include "Materials/inc/Material.hpp"

namespace materials
{
	class Heat : public Material
	{
	public:
		//constructors
		Heat(void);

		//destructor
		~Heat(void) override;

		//serialization
		void load(FILE*) override;
		void save(FILE*) const override;

		//data
		double capacity(double);
		double capacity(void) const;
		
		double condutivity(double);
		double condutivity(void) const;

	protected:
		//data
		double m_capacity;
		double m_condutivity;

		//friends
		friend class Material;
	};
}