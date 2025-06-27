#pragma once

//materials
#include "Materials/Materials/inc/Mechanic/Mechanic.hpp"

namespace materials
{
	class Steel : public Mechanic
	{
	protected:
		//constructors
		Steel(void);

		//destructor
		~Steel(void) override;

		//serialization
		void load(FILE*) override;
		void save(FILE*) const override;

	public:
		//type
		materials::type type(void) const;

		//data
		double yield_stress(double);
		double yield_stress(void) const;

		double plastic_modulus(double);
		double plastic_modulus(void) const;

		double poisson_ratio(double);
		double poisson_ratio(void) const override;

		double elastic_modulus(double);
		double elastic_modulus(void) const override;

		double thermal_coefficient(void) const override;

		//plastic stress
		double plastic_stress(double) const override;
		double plastic_gradient(double) const override;

		//yield criteria
		double yield_function(const double*, uint32_t) const override;
		void yield_hessian(double*, const double*, uint32_t) const override;
		void yield_gradient(double*, const double*, uint32_t) const override;

	protected:
		//data
		double m_yield_stress;
		double m_poisson_ratio;
		double m_plastic_modulus;
		double m_elastic_modulus;

		//friends
		friend class Material;
	};
}