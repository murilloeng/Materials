#pragma once

//Materials
#include "Materials/inc/Mechanic/Multiaxial.hpp"

namespace materials
{
	class Concrete : public Multiaxial
	{
	public:
		//type
		enum class Model : uint32_t
		{
			drucker_prager,
			bresler_pister
		};

		//constructors
		Concrete(void);

		//destructor
		~Concrete(void) override;

		//serialization
		void load(FILE*) override;
		void save(FILE*) const override;

		//type
		Type type(void) const override;

		//data
		Model model(Model);
		Model model(void) const;

		double poisson_ratio(void) const override;

		double elastic_modulus(double);
		double elastic_modulus(void) const override;

		double yield_stress_tension(double);
		double yield_stress_tension(void) const;

		double yield_stress_biaxial(double);
		double yield_stress_biaxial(void) const;

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
		Model m_model;
		double m_elastic_modulus;
		double m_yield_stress_tension;
		double m_yield_stress_biaxial;
		double m_yield_stress_compression;
	};
}