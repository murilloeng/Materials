//std
#include <cmath>
#include <cstring>

//Math
#include "Math/inc/Miscellaneous/bits.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

//Materials
#include "Materials/inc/Mechanic/Stress.hpp"

namespace materials
{
	//constructor
	Stress::Stress(const double* data, uint32_t types) : m_types{types}, m_data{data}
	{
		return;
	}

	//destructor
	Stress::~Stress(void)
	{
		return;
	}

	//operators
	double Stress::operator()(Type type) const
	{
		const uint32_t utype = 1 << uint32_t(type);
		return m_types & utype ? m_data[math::bit_index(m_types, utype)] : 0;
	}
	double Stress::operator[](Type type) const
	{
		const uint32_t utype = 1 << uint32_t(type);
		return m_types & utype ? m_data[math::bit_index(m_types, utype)] : 0;
	}

	//invariants
	double Stress::invariant_I1(void) const
	{
		//stress
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		//return
		return s11 + s22 + s33;
	}
	double Stress::invariant_I2(void) const
	{
		//stress
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//return
		return s11 * s22 + s11 * s33 + s22 * s33 - s12 * s12 - s13 * s13 - s23 * s23;
	}
	double Stress::invariant_I3(void) const
	{
		//stress
		const double s11 = (*this)[Type ::s11];
		const double s22 = (*this)[Type ::s22];
		const double s33 = (*this)[Type ::s33];
		const double s12 = (*this)[Type ::s12];
		const double s13 = (*this)[Type ::s13];
		const double s23 = (*this)[Type ::s23];
		//return
		return s11 * s22 * s33 + 2 * s12 * s13 * s23 - s11 * s23 * s23 - s22 * s13 * s13 - s33 * s12 * s12;
	}

	//principal
	double Stress::principal_s1(void) const
	{
		//data
		const double tl = lode_angle();
		const double sm = hydrostatic_stress();
		const double J2 = deviatoric_invariant_J2();
		//return
		return sm + 2 * sqrt(J2 / 3) * cos(tl + 2 * M_PI / 3);
	}
	double Stress::principal_s2(void) const
	{
		//data
		const double tl = lode_angle();
		const double sm = hydrostatic_stress();
		const double J2 = deviatoric_invariant_J2();
		//return
		return sm + 2 * sqrt(J2 / 3) * cos(tl - 2 * M_PI / 3);
	}
	double Stress::principal_s3(void) const
	{
		//data
		const double tl = lode_angle();
		const double sm = hydrostatic_stress();
		const double J2 = deviatoric_invariant_J2();
		//return
		return sm + 2 * sqrt(J2 / 3) * cos(tl);
	}

	//gradients
	double* Stress::gradient_I1(double* dI1) const
	{
		//gradient
		const double dIc[] = {1, 1, 1, 0, 0, 0};
		//return
		return filter_vector(dI1, dIc);
	}
	double* Stress::gradient_I2(double* dI2) const
	{
		//data
		double dIc[6];
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//I2 gradient
		dIc[3] = -2 * s12;
		dIc[4] = -2 * s13;
		dIc[5] = -2 * s23;
		dIc[0] = s22 + s33;
		dIc[1] = s11 + s33;
		dIc[2] = s11 + s22;
		//return
		return filter_vector(dI2, dIc);
	}
	double* Stress::gradient_I3(double* dI3) const
	{
		//data
		double dIc[6];
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//I3 gradient
		dIc[0] = s22 * s33 - s23 * s23;
		dIc[1] = s11 * s33 - s13 * s13;
		dIc[2] = s11 * s22 - s12 * s12;
		dIc[3] = -2 * (s12 * s33 + s13 * s23);
		dIc[4] = -2 * (s13 * s22 + s12 * s23);
		dIc[5] = -2 * (s23 * s11 + s12 * s13);
		//return
		return filter_vector(dI3, dIc);
	}

	//hessians
	double* Stress::hessian_I2(double* ddI2) const
	{
		//data
		double ddIc[36];
		memset(ddIc, 0, sizeof(ddIc));
		//I2 hessian
		ddIc[6 * 1 + 0] = ddIc[6 * 0 + 1] = 1;
		ddIc[6 * 2 + 0] = ddIc[6 * 0 + 2] = 1;
		ddIc[6 * 2 + 1] = ddIc[6 * 1 + 2] = 1;
		ddIc[6 * 3 + 3] = ddIc[6 * 4 + 4] = ddIc[6 * 5 + 5] = -2;
		//return
		return filter_matrix(ddI2, ddIc);
	}
	double* Stress::hessian_I3(double* ddI3) const
	{
		//data
		double ddIc[36];
		memset(ddIc, 0, sizeof(ddIc));
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//I3 hessian
		ddIc[6 * 3 + 3] = -2 * s33;
		ddIc[6 * 4 + 4] = -2 * s22;
		ddIc[6 * 5 + 5] = -2 * s11;
		ddIc[6 * 1 + 0] = ddIc[6 * 0 + 1] = s33;
		ddIc[6 * 2 + 0] = ddIc[6 * 0 + 2] = s22;
		ddIc[6 * 2 + 1] = ddIc[6 * 1 + 2] = s11;
		ddIc[6 * 5 + 0] = ddIc[6 * 0 + 5] = -2 * s23;
		ddIc[6 * 4 + 1] = ddIc[6 * 1 + 4] = -2 * s13;
		ddIc[6 * 3 + 2] = ddIc[6 * 2 + 3] = -2 * s12;
		ddIc[6 * 4 + 3] = ddIc[6 * 3 + 4] = -2 * s23;
		ddIc[6 * 5 + 3] = ddIc[6 * 3 + 5] = -2 * s13;
		ddIc[6 * 5 + 4] = ddIc[6 * 4 + 5] = -2 * s12;
		//return
		return filter_matrix(ddI3, ddIc);
	}

	//hydrostatic
	double Stress::hydrostatic_stress(void) const
	{
		return invariant_I1() / 3;
	}
	double* Stress::hydrostatic_gradient(double* dIm) const
	{
		//gradient
		const double dIc[] = {1, 1, 1, 0, 0, 0};
		//return
		return filter_vector(dIm, dIc, 3);
	}

	//deviatoric
	double Stress::deviatoric_invariant_J2(void) const
	{
		//data
		const double I1 = invariant_I1();
		const double I2 = invariant_I2();
		//return
		return I1 * I1 / 3 - I2;
	}
	double Stress::deviatoric_invariant_J3(void) const
	{
		//data
		const double I1 = invariant_I1();
		const double I2 = invariant_I2();
		const double I3 = invariant_I3();
		//return
		return 2 * I1 * I1 * I1 / 27 - I1 * I2 / 3 + I3;
	}

	double Stress::deviatoric_principal_s1(void) const
	{
		return principal_s1() - hydrostatic_stress();
	}
	double Stress::deviatoric_principal_s2(void) const
	{
		return principal_s2() - hydrostatic_stress();
	}
	double Stress::deviatoric_principal_s3(void) const
	{
		return principal_s3() - hydrostatic_stress();
	}

	double* Stress::deviatoric_gradient_J2(double* dJ2) const
	{
		//data
		double dJc[6];
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//J2 gradient
		dJc[3] = 2 * s12;
		dJc[4] = 2 * s13;
		dJc[5] = 2 * s23;
		dJc[0] = (2 * s11 - s22 - s33) / 3;
		dJc[1] = (2 * s22 - s11 - s33) / 3;
		dJc[2] = (2 * s33 - s11 - s22) / 3;
		//return
		return filter_vector(dJ2, dJc);
	}
	double* Stress::deviatoric_gradient_J3(double* dJ3) const
	{
		//data
		double dJc[6];
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//J3 gradient
		dJc[3] = 2 * s12 * (s11 + s22 - 2 * s33) / 3 + 2 * s13 * s23;
		dJc[4] = 2 * s13 * (s11 + s33 - 2 * s22) / 3 + 2 * s12 * s23;
		dJc[5] = 2 * s23 * (s22 + s33 - 2 * s11) / 3 + 2 * s12 * s13;
		dJc[0] = (2 * s11 * s11 - s22 * s22 - s33 * s33 + 4 * s22 * s33 - 2 * s11 * (s22 + s33) + 3 * (s12 * s12 + s13 * s13 - 2 * s23 * s23)) / 9;
		dJc[1] = (2 * s22 * s22 - s11 * s11 - s33 * s33 + 4 * s11 * s33 - 2 * s22 * (s11 + s33) + 3 * (s12 * s12 + s23 * s23 - 2 * s13 * s13)) / 9;
		dJc[2] = (2 * s33 * s33 - s11 * s11 - s22 * s22 + 4 * s11 * s22 - 2 * s33 * (s11 + s22) + 3 * (s13 * s13 + s23 * s23 - 2 * s12 * s12)) / 9;
		//return
		return filter_vector(dJ3, dJc);
	}

	double* Stress::deviatoric_hessian_J2(double* ddJ2) const
	{
		//data
		double ddJc[36];
		memset(ddJc, 0, sizeof(ddJc));
		//J2 hessian
		ddJc[6 * 1 + 0] = ddJc[6 * 0 + 1] = -1.0 / 3;
		ddJc[6 * 2 + 0] = ddJc[6 * 0 + 2] = -1.0 / 3;
		ddJc[6 * 2 + 1] = ddJc[6 * 1 + 2] = -1.0 / 3;
		ddJc[6 * 3 + 3] = ddJc[6 * 4 + 4] = ddJc[6 * 5 + 5] = 2;
		ddJc[6 * 0 + 0] = ddJc[6 * 1 + 1] = ddJc[6 * 2 + 2] = 2.0 / 3;
		//return
		return filter_matrix(ddJ2, ddJc);
	}
	double* Stress::deviatoric_hessian_J3(double* ddJ3) const
	{
		//data
		double ddJc[36];
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//J3 hessian
		ddJc[6 * 4 + 3] = ddJc[6 * 3 + 4] = -2 * s23;
		ddJc[6 * 5 + 3] = ddJc[6 * 3 + 5] = -2 * s13;
		ddJc[6 * 5 + 4] = ddJc[6 * 4 + 5] = -2 * s12;
		ddJc[6 * 0 + 0] = 2 * (2 * s11 - s22 - s33) / 9;
		ddJc[6 * 1 + 1] = 2 * (2 * s22 - s11 - s33) / 9;
		ddJc[6 * 2 + 2] = 2 * (2 * s33 - s11 - s22) / 9;
		ddJc[6 * 3 + 3] = 2 * (s11 + s22 - 2 * s33) / 3;
		ddJc[6 * 4 + 4] = 2 * (s11 + s33 - 2 * s22) / 3;
		ddJc[6 * 5 + 5] = 2 * (s22 + s33 - 2 * s11) / 3;
		ddJc[6 * 5 + 0] = ddJc[6 * 0 + 5] = -4 * s23 / 3;
		ddJc[6 * 4 + 1] = ddJc[6 * 1 + 4] = -4 * s13 / 3;
		ddJc[6 * 3 + 2] = ddJc[6 * 2 + 3] = -4 * s12 / 3;
		ddJc[6 * 3 + 0] = ddJc[6 * 0 + 3] = +2 * s12 / 3;
		ddJc[6 * 4 + 0] = ddJc[6 * 0 + 4] = +2 * s13 / 3;
		ddJc[6 * 3 + 1] = ddJc[6 * 1 + 3] = +2 * s12 / 3;
		ddJc[6 * 5 + 1] = ddJc[6 * 1 + 5] = +2 * s23 / 3;
		ddJc[6 * 4 + 2] = ddJc[6 * 2 + 4] = +2 * s13 / 3;
		ddJc[6 * 5 + 2] = ddJc[6 * 2 + 5] = +2 * s23 / 3;
		ddJc[6 * 1 + 0] = ddJc[6 * 0 + 1] = 2 * (2 * s33 - s11 - s22) / 9;
		ddJc[6 * 2 + 0] = ddJc[6 * 0 + 2] = 2 * (2 * s22 - s11 - s33) / 9;
		ddJc[6 * 2 + 1] = ddJc[6 * 1 + 2] = 2 * (2 * s11 - s22 - s33) / 9;
		//return
		return filter_matrix(ddJ3, ddJc);
	}

	//von mises
	double Stress::von_mises_stress(void) const
	{
		return sqrt(3 * deviatoric_invariant_J2());
	}
	double* Stress::von_mises_hessian(double* ddse) const
	{
		//data
		double ddsc[36];
		const double se = von_mises_stress();
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//von mises hessian
		ddsc[6 * 4 + 3] = ddsc[6 * 3 + 4] = -12 * s12 * s13;
		ddsc[6 * 5 + 3] = ddsc[6 * 3 + 5] = -12 * s12 * s23;
		ddsc[6 * 5 + 4] = ddsc[6 * 4 + 5] = -12 * s13 * s23;
		ddsc[6 * 3 + 0] = ddsc[6 * 0 + 3] = 2 * s12 * (s22 + s33 - 2 * s11);
		ddsc[6 * 4 + 0] = ddsc[6 * 0 + 4] = 2 * s13 * (s22 + s33 - 2 * s11);
		ddsc[6 * 5 + 0] = ddsc[6 * 0 + 5] = 2 * s23 * (s22 + s33 - 2 * s11);
		ddsc[6 * 3 + 1] = ddsc[6 * 1 + 3] = 2 * s12 * (s11 + s33 - 2 * s22);
		ddsc[6 * 4 + 1] = ddsc[6 * 1 + 4] = 2 * s13 * (s11 + s33 - 2 * s22);
		ddsc[6 * 5 + 1] = ddsc[6 * 1 + 5] = 2 * s23 * (s11 + s33 - 2 * s22);
		ddsc[6 * 3 + 2] = ddsc[6 * 2 + 3] = 2 * s12 * (s11 + s22 - 2 * s33);
		ddsc[6 * 4 + 2] = ddsc[6 * 2 + 4] = 2 * s13 * (s11 + s22 - 2 * s33);
		ddsc[6 * 5 + 2] = ddsc[6 * 2 + 5] = 2 * s23 * (s11 + s22 - 2 * s33);
		ddsc[6 * 0 + 0] = pow(s22 - s33, 2) + 4 * (s12 * s12 + s13 * s13 + s23 * s23);
		ddsc[6 * 1 + 1] = pow(s11 - s33, 2) + 4 * (s12 * s12 + s13 * s13 + s23 * s23);
		ddsc[6 * 2 + 2] = pow(s11 - s22, 2) + 4 * (s12 * s12 + s13 * s13 + s23 * s23);
		ddsc[6 * 1 + 0] = ddsc[6 * 0 + 1] = (s11 - s33) * (s33 - s22) - 2 * (s12 * s12 + s13 * s13 + s23 * s23);
		ddsc[6 * 2 + 0] = ddsc[6 * 0 + 2] = (s11 - s22) * (s22 - s33) - 2 * (s12 * s12 + s13 * s13 + s23 * s23);
		ddsc[6 * 2 + 1] = ddsc[6 * 1 + 2] = (s11 - s22) * (s33 - s11) - 2 * (s12 * s12 + s13 * s13 + s23 * s23);
		ddsc[6 * 3 + 3] = 4 * (s11 * s11 + s22 * s22 + s33 * s33 - s22 * s33 - s11 * (s22 + s33) + 3 * (s13 * s13 + s23 * s23));
		ddsc[6 * 4 + 4] = 4 * (s11 * s11 + s22 * s22 + s33 * s33 - s22 * s33 - s11 * (s22 + s33) + 3 * (s12 * s12 + s23 * s23));
		ddsc[6 * 5 + 5] = 4 * (s11 * s11 + s22 * s22 + s33 * s33 - s22 * s33 - s11 * (s22 + s33) + 3 * (s12 * s12 + s13 * s13));
		//return
		return filter_matrix(ddse, ddsc, se == 0 ? 1 : 4 * se * se * se / 3);
	}
	double* Stress::von_mises_gradient(double* dse) const
	{
		//data
		double dsc[6];
		const double se = von_mises_stress();
		const double s11 = (*this)[Type::s11];
		const double s22 = (*this)[Type::s22];
		const double s33 = (*this)[Type::s33];
		const double s12 = (*this)[Type::s12];
		const double s13 = (*this)[Type::s13];
		const double s23 = (*this)[Type::s23];
		//von misses gradient
		dsc[3] = se != 0 ? 3 * s12 / se : 1;
		dsc[4] = se != 0 ? 3 * s13 / se : 1;
		dsc[5] = se != 0 ? 3 * s23 / se : 1;
		dsc[0] = se != 0 ? (s11 - s22 / 2 - s33 / 2) / se : 1;
		dsc[1] = se != 0 ? (s22 - s11 / 2 - s33 / 2) / se : 1;
		dsc[2] = se != 0 ? (s33 - s11 / 2 - s22 / 2) / se : 1;
		//return
		return filter_vector(dse, dsc);
	}

	//angle
	double Stress::lode_angle(void) const
	{
		const double J2 = deviatoric_invariant_J2();
		const double J3 = deviatoric_invariant_J3();
		return fabs(J2) > 1e-5 ? acos(math::bound(3 * sqrt(3) / 2 * J3 / sqrt(J2 * J2 * J2))) / 3 : 0;
	}

	//filters
	double* Stress::filter_vector(double* vp, const double* vc, double s) const
	{
		for(uint32_t i = 0, c = 0; i < 6; i++)
		{
			if(m_types & 1 << i)
			{
				vp[c++] = vc[i] / s;
			}
		}
		return vp;
	}
	double* Stress::filter_matrix(double* Mp, const double* Mc, double s) const
	{
		//data
		const uint32_t n = math::bit_count(m_types);
		//filter
		for(uint32_t i = 0, a = 0; i < 6; i++)
		{
			if(m_types & 1 << i)
			{
				for(uint32_t j = 0, b = 0; j < 6; j++)
				{
					if(m_types & 1 << j)
					{
						Mp[a + n * b] = Mc[i + 6 * j] / s;
						b++;
					}
				}
				a++;
			}
		}
		return Mp;
	}
}