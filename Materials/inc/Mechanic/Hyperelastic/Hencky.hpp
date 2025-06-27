//fea
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Hyperelastic.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			class Hencky : public Hyperelastic
			{
			public:
				//constructors
				Hencky(Mechanic*);

				//destructor
				~Hencky(void);

				//type
				hyperelastic_type type(void) const override;

				//energy
				double energy(const double*) const override;
				void energy_hessian(double*, const double*) const override;
				void energy_gradient(double*, const double*) const override;
			};
		}
	}
}