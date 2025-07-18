#pragma once

//materials
#include "Materials/Materials/inc/Material.hpp"
#include "Materials/Materials/inc/Mechanic/Point.hpp"

namespace materials
{
	class Mechanic : public Material
	{
	public:
		//constructors
		Mechanic(void);

		//destructor
		virtual ~Mechanic(void);

		//serialization
		void load(FILE*) override;
		void save(FILE*) const override;

		//properties
		static bool fracture(void);
		static bool fracture(bool);

		static bool hardening(void);
		static bool hardening(bool);

		static bool inelastic(void);
		static bool inelastic(bool);
		
		double break_strain(double);
		double break_strain(void) const;

		double bulk_modulus(void) const;
		double shear_modulus(void) const;

		virtual double poisson_ratio(void) const = 0;
		virtual double elastic_modulus(void) const = 0;
		virtual double thermal_coefficient(void) const = 0;

		double* elastic_stiffness(double*, uint32_t) const;
		double* elastic_flexibility(double*, uint32_t) const;

		//plastic stress
		virtual double plastic_stress(double) const = 0;
		virtual double plastic_gradient(double) const = 0;

		//yield surface
		double yield_stress(const double*, uint32_t) const;

		//yield surface
		virtual double yield_function(const double*, uint32_t) const = 0;
		virtual void yield_hessian(double*, const double*, uint32_t) const = 0;
		virtual void yield_gradient(double*, const double*, uint32_t) const = 0;

		//return mapping
		void stress(double*, const double*, const Point&) const;
		void stiffness(double*, const double*, const Point&) const;
		void return_mapping(double*, double*, const double*, Point&) const;

		//test
		void test_plastic_gradient(void) const;
		void test_yield_hessian(uint32_t) const;
		void test_yield_gradient(uint32_t) const;

	protected:
		//return mapping
		bool broken(double*, double*, const double*, Point&) const;
		bool predictor(double*, double*, const double*, Point&) const;
		bool corrector(double*, double*, const double*, Point&) const;

		//data
		double m_break_strain;

		static bool m_break;
		static bool m_hardening;
		static bool m_inelastic;
	};
}