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
		
		double plastic_modulus(double);
		double plastic_modulus(void) const;
		
		Function yield_stress_function(Function);
		Function yield_stress_function(void) const;
		
		Function plastic_modulus_function(Function);
		Function plastic_modulus_function(void) const;

		double plastic_stress(double) const override;
		double plastic_gradient(double) const override;

		double poisson_ratio(void) const override;

		double elastic_modulus(double);
		double elastic_modulus(void) const override;

		double thermal_coefficient(double);
		double thermal_coefficient(void) const override;

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
		double m_thermal_coefficient;

		static bool m_buckling;
		std::function<double(double)> m_yield_stress_function;
		std::function<double(double)> m_plastic_modulus_function;
	};
}