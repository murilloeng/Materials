#pragma once

//materials
#include "Materials/Materials/inc/Mechanic/Mechanic.hpp"

namespace materials
{
	class Concrete : public Mechanic
	{
	public:
		//constructors
		Concrete(void);

		//destructor
		~Concrete(void) override;

		//serialization
		void load(FILE*) override;
		void save(FILE*) const override;

		//type
		materials::type type(void) const override;

		//data
		double poisson_ratio(void) const override;

		double elastic_modulus(double);
		double elastic_modulus(void) const override;

		double yield_stress_tension(double);
		double yield_stress_tension(void) const;

		double yield_stress_compression(double);
		double yield_stress_compression(void) const;

		double thermal_coefficient(void) const override;

		//plastic stress
		double plastic_stress(double) const override;
		double plastic_gradient(double) const override;

		//yield criterion
		double yield_function(const double*, uint32_t) const override;
		void yield_hessian(double*, const double*, uint32_t) const override;
		void yield_gradient(double*, const double*, uint32_t) const override;

	protected:
		//data
		double m_elastic_modulus;
		double m_yield_stress_tension;
		double m_yield_stress_compression;

		//friends
		friend class Material;
	};
}