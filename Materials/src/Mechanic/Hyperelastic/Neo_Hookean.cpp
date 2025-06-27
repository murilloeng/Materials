//std
#include <cmath>

//fea
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Mechanic.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Types.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Neo_Hookean.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			//constructors
			Neo_Hookean::Neo_Hookean(Mechanic* mechanic) : Hyperelastic(mechanic)
			{
				return;
			}

			//destructor
			Neo_Hookean::~Neo_Hookean(void)
			{
				return;
			}

			//type
			hyperelastic_type Neo_Hookean::type(void) const
			{
				return hyperelastic_type::neo_hookean;
			}

			//energy
			double Neo_Hookean::energy(const double* s) const
			{
				return 0;
			}
			void Neo_Hookean::energy_hessian(double* d2W, const double* s) const
			{
				d2W[0 + 3 * 0] = 0;
				d2W[1 + 3 * 1] = 0;
				d2W[2 + 3 * 2] = 0;
				d2W[0 + 3 * 1] = d2W[1 + 3 * 0] = 0;
				d2W[0 + 3 * 2] = d2W[2 + 3 * 0] = 0;
				d2W[1 + 3 * 2] = d2W[2 + 3 * 1] = 0;
			}
			void Neo_Hookean::energy_gradient(double* dW, const double* s) const
			{
				dW[0] = 0;
				dW[1] = 0;
				dW[2] = 0;
			}
		}
	}
}