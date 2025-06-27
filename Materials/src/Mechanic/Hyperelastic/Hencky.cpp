//std
#include <cmath>

//fea
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Mechanic.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Types.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Hencky.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			//constructors
			Hencky::Hencky(Mechanic* mechanic) : Hyperelastic(mechanic)
			{
				return;
			}

			//destructor
			Hencky::~Hencky(void)
			{
				return;
			}

			//type
			hyperelastic_type Hencky::type(void) const
			{
				return hyperelastic_type::hencky;
			}

			//energy
			double Hencky::energy(const double* s) const
			{
				//data
				const double v = m_mechanic->poisson_ratio();
				const double E = m_mechanic->elastic_modulus();
				//stretch
				const double e1 = log(s[0]);
				const double e2 = log(s[1]);
				const double e3 = log(s[2]);
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				return G * (e1 * e1 + e2 * e2 + e3 * e3) + l / 2 * am * am;
			}
			void Hencky::energy_hessian(double* d2W, const double* s) const
			{
				//data
				const double v = m_mechanic->poisson_ratio();
				const double E = m_mechanic->elastic_modulus();
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = log(s1);
				const double e2 = log(s2);
				const double e3 = log(s3);
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				d2W[0 + 3 * 1] = d2W[1 + 3 * 0] = l / s1 / s2;
				d2W[0 + 3 * 2] = d2W[2 + 3 * 0] = l / s1 / s3;
				d2W[1 + 3 * 2] = d2W[2 + 3 * 1] = l / s2 / s3;
				d2W[0 + 3 * 0] = (2 * G * (1 - e1) + l * (1 - am)) / s1 / s1;
				d2W[1 + 3 * 1] = (2 * G * (1 - e2) + l * (1 - am)) / s2 / s2;
				d2W[2 + 3 * 2] = (2 * G * (1 - e3) + l * (1 - am)) / s3 / s3;
			}
			void Hencky::energy_gradient(double* dW, const double* s) const
			{
				//data
				const double v = m_mechanic->poisson_ratio();
				const double E = m_mechanic->elastic_modulus();
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = log(s1);
				const double e2 = log(s2);
				const double e3 = log(s3);
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				dW[0] = (2 * G * e1 + l * am) / s1;
				dW[1] = (2 * G * e2 + l * am) / s2;
				dW[2] = (2 * G * e3 + l * am) / s3;
			}
		}
	}
}