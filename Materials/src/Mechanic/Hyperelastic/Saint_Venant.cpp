//std
#include <cmath>

//math
#include "Math/Math/inc/misc/bits.hpp"
#include "Math/Math/inc/linear/vector.hpp"

//fea
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/stress.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Mechanic.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Types.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Saint_Venant.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			//constructors
			Saint_Venant::Saint_Venant(Mechanic* mechanic) : Hyperelastic(mechanic)
			{
				return;
			}

			//destructor
			Saint_Venant::~Saint_Venant(void)
			{
				return;
			}

			//type
			hyperelastic_type Saint_Venant::type(void) const
			{
				return hyperelastic_type::saint_venant;
			}

			//energy
			double Saint_Venant::energy(const double* s) const
			{
				//data
				const double v = m_mechanic->poisson_ratio();
				const double E = m_mechanic->elastic_modulus();
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = (s1 * s1 - 1) / 2;
				const double e2 = (s2 * s2 - 1) / 2;
				const double e3 = (s3 * s3 - 1) / 2;
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				return G * (e1 * e1 + e2 * e2 + e3 * e3) + l / 2 * am * am;
			}
			void Saint_Venant::energy_hessian(double* d2W, const double* s) const
			{
				//data
				const double v = m_mechanic->poisson_ratio();
				const double E = m_mechanic->elastic_modulus();
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = (s1 * s1 - 1) / 2;
				const double e2 = (s2 * s2 - 1) / 2;
				const double e3 = (s3 * s3 - 1) / 2;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				d2W[0 + 3 * 1] = d2W[1 + 3 * 0] = l * s1 * s2;
				d2W[0 + 3 * 2] = d2W[2 + 3 * 0] = l * s1 * s3;
				d2W[1 + 3 * 2] = d2W[2 + 3 * 1] = l * s2 * s3;
				d2W[0 + 3 * 0] = (2 * G + l) * (3 * e1 + 1) + l * (e2 + e3);
				d2W[1 + 3 * 1] = (2 * G + l) * (3 * e2 + 1) + l * (e1 + e3);
				d2W[2 + 3 * 2] = (2 * G + l) * (3 * e3 + 1) + l * (e1 + e2);
			}
			void Saint_Venant::energy_gradient(double* dW, const double* s) const
			{
				//data
				const double v = m_mechanic->poisson_ratio();
				const double E = m_mechanic->elastic_modulus();
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = (s1 * s1 - 1) / 2;
				const double e2 = (s2 * s2 - 1) / 2;
				const double e3 = (s3 * s3 - 1) / 2;
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				dW[0] = (2 * G * e1 + l * am) * s1;
				dW[1] = (2 * G * e2 + l * am) * s2;
				dW[2] = (2 * G * e3 + l * am) * s3;
			}

			//stress
			void Saint_Venant::stress(double* sa, double* Ca, const double* ea, uint32_t ma) const
			{
				//data
				double Cf[36];
				const uint32_t na = math::bit_count(ma);
				const uint32_t nb = 6 - math::bit_count(ma);
				const uint32_t mb = (2 * uint32_t(stress::last) - 3) & ~ma;
				math::matrix Caa(na, na), Cab(na, nb), Cba(nb, na), Cbb(nb, nb); 
				//setup
				math::matrix(Cf, 6, 6).zeros();
				const double v = m_mechanic->poisson_ratio();
				const double E = m_mechanic->elastic_modulus();
				Cf[3 + 6 * 3] = Cf[4 + 6 * 4] = Cf[5 + 6 * 5] = E / (1 + v) / 2;
				Cf[0 + 6 * 1] = Cf[0 + 6 * 2] = Cf[1 + 6 * 0] = E * v / (1 - 2 * v) / (1 + v);
				Cf[1 + 6 * 2] = Cf[2 + 6 * 0] = Cf[2 + 6 * 1] = E * v / (1 - 2 * v) / (1 + v);
				Cf[0 + 6 * 0] = Cf[1 + 6 * 1] = Cf[2 + 6 * 2] = E * (1 - v) / (1 - 2 * v) / (1 + v);
				//stiffness
				stiffness_block(Caa.data(), Cf, ma, ma);
				stiffness_block(Cab.data(), Cf, ma, mb);
				stiffness_block(Cba.data(), Cf, mb, ma);
				stiffness_block(Cbb.data(), Cf, mb, mb);
				math::matrix(Ca, na, na) = Caa - Cab * Cbb.inverse() * Cba;
				//stress
				math::vector(sa, na) = math::matrix(Ca, na, na) * math::vector(ea, na);
			}
		}
	}
}