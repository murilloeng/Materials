#pragma once

//std
#include <cstdint>

//materials
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Mechanic.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			enum class model : uint32_t
			{
				gent,
				hencky,
				neo_hookean,
				saint_venant
			};

			class Hyperelastic2 : public Mechanic
			{
			public:
				//friends
				friend class Material;

			protected:
				//constructors
				Hyperelastic2(void);

				//destructor
				~Hyperelastic2(void);

				//serialization
				void load(FILE*) override;
				void save(FILE*) const override;

			public:
				//name
				const char* model_name(void) const;
				static const char* model_name(model);

				//type
				materials::type type(void) const;

				//data
				double poisson_ratio(double);
				double poisson_ratio(void) const override;

				double elastic_modulus(double);
				double elastic_modulus(void) const override;

				double thermal_coefficient(double);
				double thermal_coefficient(void) const override;

				materials::model model(void) const;
				materials::model model(materials::model);

				//plastic stress
				double plastic_stress(double) const override;
				double plastic_gradient(double) const override;

				//yield criteria
				double yield_function(const double*, uint32_t) const override;
				void yield_hessian(double*, const double*, uint32_t) const override;
				void yield_gradient(double*, const double*, uint32_t) const override;

				//energy
				double energy(const double*) const;
				void energy_hessian(double*, const double*) const;
				void energy_gradient(double*, const double*) const;

				//states
				void plane_stress(double*, double*, double*) const;
				void uniaxial_stress(double*, double*, double*) const;
			
			private:
				//energy
				double energy_gent(const double*) const;
				void energy_gent_hessian(double*, const double*) const;
				void energy_gent_gradient(double*, const double*) const;

				double energy_hencky(const double*) const;
				void energy_hencky_hessian(double*, const double*) const;
				void energy_hencky_gradient(double*, const double*) const;

				double energy_neo_hookean(const double*) const;
				void energy_neo_hookean_hessian(double*, const double*) const;
				void energy_neo_hookean_gradient(double*, const double*) const;

				double energy_saint_venant(const double*) const;
				void energy_saint_venant_hessian(double*, const double*) const;
				void energy_saint_venant_gradient(double*, const double*) const;

				//data
				materials::model m_model;

				double m_poisson_ratio;
				double m_elastic_modulus;
				double m_thermal_coefficient;
			};
		}
	}
}