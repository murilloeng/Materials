#pragma once

//std
#include <functional>

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

		//types
		typedef std::function<double(double)> Function;

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

		Function yield_stress_function(Function);
		Function yield_stress_function(void) const;

		Function plastic_modulus_function(Function);
		Function plastic_modulus_function(void) const;

		//return mapping
		void stress(double&, double, const Point&) const;
		void stiffness(double&, double, const Point&) const;
		void return_mapping(double&, double&, double, Point&) const;

	protected:
		//compute
		double compute_yield_stress(double) const;
		double compute_plastic_modulus(double) const;

		//data
		double m_yield_stress;
		double m_buckling_stress;
		double m_elastic_modulus;
		double m_plastic_modulus;

		static bool m_buckling;
		std::function<double(double)> m_yield_stress_function;
		std::function<double(double)> m_plastic_modulus_function;
	};
}