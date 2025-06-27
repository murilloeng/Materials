//std
#include <cmath>

//fea
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Mechanic.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Gent.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Types.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			//constructors
			Gent::Gent(Mechanic* mechanic) : Hyperelastic(mechanic)
			{
				return;
			}

			//destructor
			Gent::~Gent(void)
			{
				return;
			}

			//type
			hyperelastic_type Gent::type(void) const
			{
				return hyperelastic_type::gent;
			}

			//energy
			double Gent::energy(const double* s) const
			{
				return 0;
			}
			void Gent::energy_hessian(double* d2W, const double* s) const
			{
				d2W[0 + 3 * 0] = 0;
				d2W[1 + 3 * 1] = 0;
				d2W[2 + 3 * 2] = 0;
				d2W[0 + 3 * 1] = d2W[1 + 3 * 0] = 0;
				d2W[0 + 3 * 2] = d2W[2 + 3 * 0] = 0;
				d2W[1 + 3 * 2] = d2W[2 + 3 * 1] = 0;
			}
			void Gent::energy_gradient(double* dW, const double* s) const
			{
				dW[0] = 0;
				dW[1] = 0;
				dW[2] = 0;
			}
		}
	}
}