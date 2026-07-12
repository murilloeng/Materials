#pragma once

//Materials
#include "Materials/inc/Mechanic/Mechanic.hpp"

namespace materials
{
	class Uniaxial : public Mechanic
	{
	public:
		//constructors
		Uniaxial(void);

		//destructor
		virtual ~Uniaxial(void);

		//data
		static bool buckling(void);
		static bool buckling(bool);

		double yield_stress(double);
		double yield_stress(void) const;

		double buckling_stress(double);
		double buckling_stress(void) const;

		double elastic_modulus(double);
		double elastic_modulus(void) const;

		double plastic_modulus(double);
		double plastic_modulus(void) const;

		//return mapping
		void stress(double&, double, const Point&) const;
		void stiffness(double&, double, const Point&) const;
		void return_mapping(double&, double&, double, Point&) const;

	protected:
		//data
		double m_yield_stress;
		double m_buckling_stress;
		double m_elastic_modulus;
		double m_plastic_modulus;

		static bool m_buckling;
	};
}