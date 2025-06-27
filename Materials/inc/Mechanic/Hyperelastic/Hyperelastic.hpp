#pragma once

//std
#include <cstdint>

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			class Mechanic;
			enum class hyperelastic_type : uint32_t;
		}
	}
}

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			class Hyperelastic
			{
			public:
				//constructors
				Hyperelastic(const Mechanic*);

				//destructor
				~Hyperelastic(void);

				//type
				const char* name(void) const;
				static const char* name(hyperelastic_type);
				virtual hyperelastic_type type(void) const = 0;

				//energy
				virtual double energy(const double*) const = 0;
				virtual void energy_hessian(double*, const double*) const = 0;
				virtual void energy_gradient(double*, const double*) const = 0;

				//stress
				virtual void stress(double*, double*, const double*, uint32_t) const;

			protected:
				//energy
				void quadratic_hessian(double*, const double*) const;
				void quadratic_gradient(double*, const double*) const;

				//stress
				void stress_full(double*, const double*, const double*) const;
				void stiffness_full(double*, const double*, const double*) const;

				static void stress_block(double*, const double*, uint32_t);
				static void stress_entry(double*, const double*, uint32_t);

				static void stiffness(double*, const double*, uint32_t);
				static void stiffness_block(double*, const double*, uint32_t, uint32_t);

				//kinematics
				static void kinematics(double*, double*, const double*);

				//data
				const Mechanic* m_mechanic;
			};
		}
	}
}