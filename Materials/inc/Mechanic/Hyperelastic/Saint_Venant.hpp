//fea
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Hyperelastic.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			class Saint_Venant : public Hyperelastic
			{
			public:
				//constructors
				Saint_Venant(Mechanic*);

				//destructor
				~Saint_Venant(void);

				//type
				hyperelastic_type type(void) const override;

				//energy
				double energy(const double*) const override;
				void energy_hessian(double*, const double*) const override;
				void energy_gradient(double*, const double*) const override;

				//stress
				void stress(double*, double*, const double*, uint32_t) const override;
			};
		}
	}
}