#pragma once

//Materials
#include "Materials/inc/Mechanic/Mechanic.hpp"

namespace materials
{
	class Multiaxial : public Mechanic
	{
	public:
		//constructors
		Multiaxial(void);

		//destructor
		virtual ~Multiaxial(void);

		//elasticity
		double* elastic_stiffness(double*, uint32_t) const;
		double* elastic_flexibility(double*, uint32_t) const;

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
	};
}