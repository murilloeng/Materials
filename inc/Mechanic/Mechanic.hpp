#pragma once

//Materials
#include "Materials/inc/Material.hpp"
#include "Materials/inc/Mechanic/Point.hpp"

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

		//data
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

		virtual double plastic_stress(double) const = 0;
		virtual double plastic_gradient(double) const = 0;

		virtual double poisson_ratio(void) const = 0;
		virtual double elastic_modulus(void) const = 0;
		virtual double thermal_coefficient(void) const = 0;

	protected:
		//data
		double m_break_strain;

		static bool m_break;
		static bool m_hardening;
		static bool m_inelastic;
	};
}