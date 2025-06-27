//std
#include <cmath>
#include <malloc.h>

//math
#include "Math/Math/inc/misc/bits.hpp"
#include "Math/Math/inc/linear/vector.hpp"

//fea
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/stress.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Types.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic/Hyperelastic.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			//constructors
			Hyperelastic::Hyperelastic(const Mechanic* mechanic) : m_mechanic(mechanic)
			{
				return;
			}

			//destructor
			Hyperelastic::~Hyperelastic(void)
			{
				return;
			}

			//type
			const char* Hyperelastic::name(void) const
			{
				return name(type());
			}
			const char* Hyperelastic::name(hyperelastic_type type)
			{
				hyperelastic_type types[] = {
					hyperelastic_type::gent, 
					hyperelastic_type::hencky, 
					hyperelastic_type::neo_hookean,
					hyperelastic_type::saint_venant
				};
				const char* names[] = {"Gent", "Hencky", "Neo Hookean", "Saint Venant"};
				for(uint32_t i = 0; i < sizeof(names) / sizeof(const char*); i++)
				{
					if(type == types[i]) return names[i];
				}
				return "Error";
			}

			//stress
			void Hyperelastic::stress(double* sa, double* Ca, const double* ea, uint32_t ma) const
			{
				//data
				const uint32_t na = math::bit_count(ma);
				const uint32_t nb = 6 - math::bit_count(ma);
				const uint32_t mb = (2 * uint32_t(stress::last) - 3) & ~ma;
				double s[3], N[9], eb[6], sb[6], Cbb[36], ef[6], sf[6], Cf[36];
				//setup
				stress_entry(ef, ea, ma);
				math::vector(eb, nb).zeros();
				math::vector ebm(eb, nb), sbm(sb, nb);
				//stress
				while(true)
				{
					//strain
					stress_entry(ef, eb, mb);
					//kinematics
					kinematics(s, N, ef);
					//stress
					stress_full(sf, s, N);
					stress_block(sa, sf, ma);
					stress_block(sb, sf, mb);
					//check
					stiffness_full(Cf, s, N);
					if(math::vector(sb, nb).norm() <= 1e-5 * math::vector(sa, na).norm())
					{
						stiffness(Ca, Cf, ma);
						break;
					}
					//update
					stiffness_block(Cbb, Cf, mb, mb);
					math::matrix(Cbb, nb, nb).solve(sbm, sbm), ebm -= sbm;
				}
			}

			//energy
			void Hyperelastic::quadratic_hessian(double* d2Wq, const double* s) const
			{
				double dW[3];
				energy_gradient(dW, s);
				energy_hessian(d2Wq, s);
				for(uint32_t k = 0; k < 3; k++)
				{
					for(uint32_t p = 0; p < 3; p++)
					{
						d2Wq[k + 3 * p] /= 4 * s[k] * s[p];
					}
					d2Wq[k + 3 * k] -= dW[k] / s[k] / s[k] / s[k] / 4;
				}
			}
			void Hyperelastic::quadratic_gradient(double* dWq, const double* s) const
			{
				energy_gradient(dWq, s);
				for(uint32_t k = 0; k < 3; k++)
				{
					dWq[k] /= 2 * s[k];
				}
			}

			//stress
			void Hyperelastic::stress_full(double* sv, const double* s, const double* N) const
			{
				double dWq[3];
				quadratic_gradient(dWq, s);
				math::vector(sv, 6).zeros();
				for(uint32_t k = 0; k < 3; k++)
				{
					const double* Nk = N + 3 * k;
					sv[0] += 2 * dWq[k] * Nk[0] * Nk[0];
					sv[1] += 2 * dWq[k] * Nk[1] * Nk[1];
					sv[2] += 2 * dWq[k] * Nk[2] * Nk[2];
					sv[3] += 2 * dWq[k] * Nk[0] * Nk[1];
					sv[4] += 2 * dWq[k] * Nk[0] * Nk[2];
					sv[5] += 2 * dWq[k] * Nk[1] * Nk[2];
				}
			}
			void Hyperelastic::stiffness_full(double* C, const double* s, const double* N) const
			{
				double km, kg;
				double dWq[3], d2Wq[9];
				quadratic_hessian(d2Wq, s);
				quadratic_gradient(dWq, s);
				math::matrix(C, 6, 6).zeros();
				const uint32_t p1[] = {0, 1, 2, 0, 0, 1};
				const uint32_t p2[] = {0, 1, 2, 1, 2, 2};
				for(uint32_t i = 0; i < 6; i++)
				{
					for(uint32_t j = 0; j < 6; j++)
					{
						for(uint32_t k = 0; k < 3; k++)
						{
							for(uint32_t p = 0; p < 3; p++)
							{
								km = 4 * d2Wq[k + 3 * p];
								const double* Nk = N + 3 * k;
								const double* Np = N + 3 * p;
								C[i + 6 * j] += km * Nk[p1[i]] * Nk[p2[i]] * Np[p1[j]] * Np[p2[j]];
								if(k != p)
								{
									if(fabs(s[k] - s[p]) < 1e-5)
									{
										kg = 4 * (d2Wq[k + 3 * k] - d2Wq[k + 3 * p]);
									}
									else
									{
										kg = 4 * (dWq[k] - dWq[p]) / (s[k] * s[k] - s[p] * s[p]);
									}
									C[i + 6 * j] += kg * Nk[p1[i]] * Np[p2[i]] * Nk[p1[j]] * Np[p2[j]];
								}
							}
						}
						if(j > 2) C[i + 6 * j] /= 2;
					}
				}
			}

			void Hyperelastic::stress_block(double* s, const double* sf, uint32_t mr)
			{
				for(uint32_t i = 0, a = 0; i < 6; i++)
				{
					if(mr & 1 << i)
					{
						s[a] = sf[i];
						a++;
					}
				}
			}
			void Hyperelastic::stress_entry(double* sf, const double* s, uint32_t mr)
			{
				for(uint32_t i = 0, a = 0; i < 6; i++)
				{
					if(mr & 1 << i)
					{
						sf[i] = s[a];
						a++;
					}
				}
			}

			void Hyperelastic::stiffness(double* Ca, const double* Cf, uint32_t ma)
			{
				//data
				const uint32_t na = math::bit_count(ma);
				const uint32_t mb = (2 * uint32_t(stress::last) - 3) & ~ma;
				math::matrix Caa(na, na), Cab(na, 6 - na), Cba(6 - na, na), Cbb(6 - na, 6 - na);
				//stiffness
				if(ma == 2 * uint32_t(stress::last) - 3)
				{
					math::matrix(Ca, na, na) = math::matrix(Cf, na, na);
				}
				else
				{
					stiffness_block(Caa.data(), Cf, ma, ma);
					stiffness_block(Cab.data(), Cf, ma, mb);
					stiffness_block(Cba.data(), Cf, mb, ma);
					stiffness_block(Cbb.data(), Cf, mb, mb);
					math::matrix(Ca, na, na) = Caa - Cab * Cbb.inverse() * Cba;
				}
			}
			void Hyperelastic::stiffness_block(double* C, const double* Cf, uint32_t mr, uint32_t mc)
			{
				const uint32_t na = math::bit_count(mr);
				for(uint32_t i = 0, a = 0; i < 6; i++)
				{
					if(mr & 1 << i)
					{
						for(uint32_t j = 0, b = 0; j < 6; j++)
						{
							if(mc & 1 << j)
							{
								C[a + na * b] = Cf[i + 6 * j];
								b++;
							}
						}
						a++;
					}
				}
			}

			//kinematics
			void Hyperelastic::kinematics(double* s, double* N, const double* ef)
			{
				//data
				math::vector sm(s, 3);
				math::matrix C(3, 3), Nm(N, 3, 3);
				//strain
				C[0 + 3 * 0] = 2 * ef[0] + 1;
				C[1 + 3 * 1] = 2 * ef[1] + 1;
				C[2 + 3 * 2] = 2 * ef[2] + 1;
				C[0 + 3 * 1] = C[1 + 3 * 0] = ef[3];
				C[0 + 3 * 2] = C[2 + 3 * 0] = ef[4];
				C[1 + 3 * 2] = C[2 + 3 * 1] = ef[5];
				//eigen
				C.eigen_sym(sm, Nm);
				sm[0] = sqrt(sm[0]);
				sm[1] = sqrt(sm[1]);
				sm[2] = sqrt(sm[2]);
			}
		}
	}
}