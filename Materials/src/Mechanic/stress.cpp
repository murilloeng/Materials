//std
#include <cmath>

//math
#include "Math/Math/inc/misc/bits.hpp"
#include "Math/Math/inc/misc/util.hpp"

//materials
#include "Materials/Materials/inc/Mechanic/stress.hpp"

namespace materials
{

	double stress_value(const double* s, uint32_t t, uint32_t c)
	{
		return t & c ? s[math::bit_index(t, c)] : 0;
	}
	double* stress_vector(double* vp, const double* vc, uint32_t t, double p)
	{
		for(uint32_t i = 0, c = 0; i < 6; i++)
		{
			if(t & 1 << i)
			{
				vp[c++] = vc[i] / p;
			}
		}
		return vp;
	}
	double* stress_matrix(double* mp, const double* mc, uint32_t t, double p)
	{
		const uint32_t n = math::bit_count(t);
		for(uint32_t i = 0, a = 0; i < 6; i++)
		{
			if(t & 1 << i)
			{
				for(uint32_t j = 0, b = 0; j < 6; j++)
				{
					if(t & 1 << j)
					{
						mp[a + n * b] = mc[i + 6 * j] / p;
						b++;
					}
				}
				a++;
			}
		}
		return mp;
	}

	double stress_invariant_I1(const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		//return
		return s11 + s22 + s33;
	}
	double stress_invariant_I2(const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//return
		return s11 * s22 + s11 * s33 + s22 * s33 - s12 * s12 - s13 * s13 - s23 * s23;
	}
	double stress_invariant_I3(const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//return
		return s11 * s22 * s33 + 2 * s12 * s13 * s23 - s11 * s23 * s23 - s22 * s13 * s13 - s33 * s12 * s12;
	}

	double stress_principal_s1(const double* s, uint32_t t)
	{
		const double tl = lode_angle(s, t);
		const double sm = hydrostatic_stress(s, t);
		const double J2 = deviatoric_invariant_J2(s, t);
		return sm + 2 * sqrt(J2 / 3) * cos(tl + 2 * M_PI / 3);
	}
	double stress_principal_s2(const double* s, uint32_t t)
	{
		const double tl = lode_angle(s, t);
		const double sm = hydrostatic_stress(s, t);
		const double J2 = deviatoric_invariant_J2(s, t);
		return sm + 2 * sqrt(J2 / 3) * cos(tl - 2 * M_PI / 3);
	}
	double stress_principal_s3(const double* s, uint32_t t)
	{
		const double tl = lode_angle(s, t);
		const double sm = hydrostatic_stress(s, t);
		const double J2 = deviatoric_invariant_J2(s, t);
		return sm + 2 * sqrt(J2 / 3) * cos(tl);
	}

	double* stress_gradient_I1(double* dI1, const double*, uint32_t t)
	{
		//I1 gradient
		const double dIc[] = {1, 1, 1, 0, 0, 0};
		//return
		return stress_vector(dI1, dIc, t);
	}
	double* stress_gradient_I2(double* dI2, const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//I2 gradient
		double dIc[6];
		dIc[3] = -2 * s12;
		dIc[4] = -2 * s13;
		dIc[5] = -2 * s23;
		dIc[0] = s22 + s33;
		dIc[1] = s11 + s33;
		dIc[2] = s11 + s22;
		//return
		return stress_vector(dI2, dIc, t);
	}
	double* stress_gradient_I3(double* dI3, const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//I3 gradient
		double dIc[6];
		dIc[0] = s22 * s33 - s23 * s23;
		dIc[1] = s11 * s33 - s13 * s13;
		dIc[2] = s11 * s22 - s12 * s12;
		dIc[3] = -2 * (s12 * s33 + s13 * s23);
		dIc[4] = -2 * (s13 * s22 + s12 * s23);
		dIc[5] = -2 * (s23 * s11 + s12 * s13);
		//return
		return stress_vector(dI3, dIc, t);
	}

	double* stress_hessian_I2(double* d2I2, const double*, uint32_t t)
	{
		//I2 hessian
		double d2Ic[36];
		d2Ic[6 * 0 + 0] = 0;
		d2Ic[6 * 1 + 1] = 0;
		d2Ic[6 * 2 + 2] = 0;
		d2Ic[6 * 3 + 3] = -2;
		d2Ic[6 * 4 + 4] = -2;
		d2Ic[6 * 5 + 5] = -2;
		d2Ic[6 * 1 + 0] = d2Ic[6 * 0 + 1] = 1;
		d2Ic[6 * 2 + 0] = d2Ic[6 * 0 + 2] = 1;
		d2Ic[6 * 2 + 1] = d2Ic[6 * 1 + 2] = 1;
		d2Ic[6 * 3 + 0] = d2Ic[6 * 0 + 3] = 0;
		d2Ic[6 * 4 + 0] = d2Ic[6 * 0 + 4] = 0;
		d2Ic[6 * 5 + 0] = d2Ic[6 * 0 + 5] = 0;
		d2Ic[6 * 3 + 1] = d2Ic[6 * 1 + 3] = 0;
		d2Ic[6 * 4 + 1] = d2Ic[6 * 1 + 4] = 0;
		d2Ic[6 * 5 + 1] = d2Ic[6 * 1 + 5] = 0;
		d2Ic[6 * 3 + 2] = d2Ic[6 * 2 + 3] = 0;
		d2Ic[6 * 4 + 2] = d2Ic[6 * 2 + 4] = 0;
		d2Ic[6 * 5 + 2] = d2Ic[6 * 2 + 5] = 0;
		d2Ic[6 * 4 + 3] = d2Ic[6 * 3 + 4] = 0;
		d2Ic[6 * 5 + 3] = d2Ic[6 * 3 + 5] = 0;
		d2Ic[6 * 5 + 4] = d2Ic[6 * 4 + 5] = 0;
		//return
		return stress_matrix(d2I2, d2Ic, t);
	}
	double* stress_hessian_I3(double* d2I3, const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//I3 hessian
		double d2Ic[36];
		d2Ic[6 * 0 + 0] = 0;
		d2Ic[6 * 1 + 1] = 0;
		d2Ic[6 * 2 + 2] = 0;
		d2Ic[6 * 3 + 3] = -2 * s33;
		d2Ic[6 * 4 + 4] = -2 * s22;
		d2Ic[6 * 5 + 5] = -2 * s11;
		d2Ic[6 * 3 + 0] = d2Ic[6 * 0 + 3] = 0;
		d2Ic[6 * 4 + 0] = d2Ic[6 * 0 + 4] = 0;
		d2Ic[6 * 3 + 1] = d2Ic[6 * 1 + 3] = 0;
		d2Ic[6 * 5 + 1] = d2Ic[6 * 1 + 5] = 0;
		d2Ic[6 * 4 + 2] = d2Ic[6 * 2 + 4] = 0;
		d2Ic[6 * 5 + 2] = d2Ic[6 * 2 + 5] = 0;
		d2Ic[6 * 1 + 0] = d2Ic[6 * 0 + 1] = s33;
		d2Ic[6 * 2 + 0] = d2Ic[6 * 0 + 2] = s22;
		d2Ic[6 * 2 + 1] = d2Ic[6 * 1 + 2] = s11;
		d2Ic[6 * 5 + 0] = d2Ic[6 * 0 + 5] = -2 * s23;
		d2Ic[6 * 4 + 1] = d2Ic[6 * 1 + 4] = -2 * s13;
		d2Ic[6 * 3 + 2] = d2Ic[6 * 2 + 3] = -2 * s12;
		d2Ic[6 * 4 + 3] = d2Ic[6 * 3 + 4] = -2 * s23;
		d2Ic[6 * 5 + 3] = d2Ic[6 * 3 + 5] = -2 * s13;
		d2Ic[6 * 5 + 4] = d2Ic[6 * 4 + 5] = -2 * s12;
		//return
		return stress_matrix(d2I3, d2Ic, t);
	}

	double hydrostatic_stress(const double* s, uint32_t t)
	{
		return stress_invariant_I1(s, t) / 3;
	}
	double* hydrostatic_gradient(double* dsm, const double*, uint32_t t)
	{
		const double dIc[] = {1, 1, 1, 0, 0, 0};
		return stress_vector(dsm, dIc, t, 3);
	}

	double deviatoric_invariant_J2(const double* s, uint32_t t)
	{
		const double I1 = stress_invariant_I1(s, t);
		const double I2 = stress_invariant_I2(s, t);
		return I1 * I1 / 3 - I2;
	}
	double deviatoric_invariant_J3(const double* s, uint32_t t)
	{
		const double I1 = stress_invariant_I1(s, t);
		const double I2 = stress_invariant_I2(s, t);
		const double I3 = stress_invariant_I3(s, t);
		return 2 * I1 * I1 * I1 / 27 - I1 * I2 / 3 + I3;
	}
	double deviatoric_principal_s1(const double* s, uint32_t t)
	{
		return stress_principal_s1(s, t) - hydrostatic_stress(s, t);
	}
	double deviatoric_principal_s2(const double* s, uint32_t t)
	{
		return stress_principal_s2(s, t) - hydrostatic_stress(s, t);
	}
	double deviatoric_principal_s3(const double* s, uint32_t t)
	{
		return stress_principal_s3(s, t) - hydrostatic_stress(s, t);
	}

	double* deviatoric_gradient_J2(double* dJ2, const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//J2 gradient
		double dJc[6];
		dJc[3] = 2 * s12;
		dJc[4] = 2 * s13;
		dJc[5] = 2 * s23;
		dJc[0] = (2 * s11 - s22 - s33) / 3;
		dJc[1] = (2 * s22 - s11 - s33) / 3;
		dJc[2] = (2 * s33 - s11 - s22) / 3;
		//return
		return stress_vector(dJ2, dJc, t);
	}
	double* deviatoric_gradient_J3(double* dJ3, const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//J3 gradient
		double dJc[6];
		dJc[3] = (2 * s12 * (s11 + s22 - 2 * s33)) / 3 - 2 * s13 * s23;
		dJc[4] = (2 * s13 * (s11 + s33 - 2 * s22)) / 3 - 2 * s12 * s23;
		dJc[5] = (2 * s23 * (s22 + s33 - 2 * s11)) / 3 - 2 * s12 * s13;
		dJc[0] = (2 * s11 * s11 - s22 * s22 - s33 * s33 + 4 * s22 * s33 - 2 * s11 * (s22 + s33) + 3 * (s12 * s12 + s13 * s13 - 2 * s23 * s23)) / 9;
		dJc[1] = (2 * s22 * s22 - s11 * s11 - s33 * s33 + 4 * s11 * s33 - 2 * s22 * (s11 + s33) + 3 * (s12 * s12 + s23 * s23 - 2 * s13 * s13)) / 9;
		dJc[2] = (2 * s33 * s33 - s11 * s11 - s22 * s22 + 4 * s11 * s22 - 2 * s33 * (s11 + s22) + 3 * (s13 * s13 + s23 * s23 - 2 * s12 * s12)) / 9;
		//return
		return stress_vector(dJ3, dJc, t);
	}

	double* deviatoric_hessian_J2(double* d2J2, const double*, uint32_t t)
	{
		//J2 hessian
		double d2Jc[36];
		d2Jc[6 * 3 + 3] = 2;
		d2Jc[6 * 4 + 4] = 2;
		d2Jc[6 * 5 + 5] = 2;
		d2Jc[6 * 0 + 0] = (double) 2 / 3;
		d2Jc[6 * 1 + 1] = (double) 2 / 3;
		d2Jc[6 * 2 + 2] = (double) 2 / 3;
		d2Jc[6 * 3 + 0] = d2Jc[6 * 0 + 3] = 0;
		d2Jc[6 * 4 + 0] = d2Jc[6 * 0 + 4] = 0;
		d2Jc[6 * 5 + 0] = d2Jc[6 * 0 + 5] = 0;
		d2Jc[6 * 3 + 1] = d2Jc[6 * 1 + 3] = 0;
		d2Jc[6 * 4 + 1] = d2Jc[6 * 1 + 4] = 0;
		d2Jc[6 * 5 + 1] = d2Jc[6 * 1 + 5] = 0;
		d2Jc[6 * 3 + 2] = d2Jc[6 * 2 + 3] = 0;
		d2Jc[6 * 4 + 2] = d2Jc[6 * 2 + 4] = 0;
		d2Jc[6 * 5 + 2] = d2Jc[6 * 2 + 5] = 0;
		d2Jc[6 * 4 + 3] = d2Jc[6 * 3 + 4] = 0;
		d2Jc[6 * 5 + 3] = d2Jc[6 * 3 + 5] = 0;
		d2Jc[6 * 5 + 4] = d2Jc[6 * 4 + 5] = 0;
		d2Jc[6 * 1 + 0] = d2Jc[6 * 0 + 1] = -(double) 1 / 3;
		d2Jc[6 * 2 + 0] = d2Jc[6 * 0 + 2] = -(double) 1 / 3;
		d2Jc[6 * 2 + 1] = d2Jc[6 * 1 + 2] = -(double) 1 / 3;
		//return
		return stress_matrix(d2J2, d2Jc, t);
	}
	double* deviatoric_hessian_J3(double* d2J3, const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//J3 hessian
		double d2Jc[36];
		d2Jc[6 * 4 + 3] = d2Jc[6 * 3 + 4] = -2 * s23;
		d2Jc[6 * 5 + 3] = d2Jc[6 * 3 + 5] = -2 * s13;
		d2Jc[6 * 5 + 4] = d2Jc[6 * 4 + 5] = -2 * s12;
		d2Jc[6 * 0 + 0] = 2 * (2 * s11 - s22 - s33) / 9;
		d2Jc[6 * 1 + 1] = 2 * (2 * s22 - s11 - s33) / 9;
		d2Jc[6 * 2 + 2] = 2 * (2 * s33 - s11 - s22) / 9;
		d2Jc[6 * 3 + 3] = 2 * (s11 + s22 - 2 * s33) / 3;
		d2Jc[6 * 4 + 4] = 2 * (s11 + s33 - 2 * s22) / 3;
		d2Jc[6 * 5 + 5] = 2 * (s22 + s33 - 2 * s11) / 3;
		d2Jc[6 * 5 + 0] = d2Jc[6 * 0 + 5] = -4 * s23 / 3;
		d2Jc[6 * 4 + 1] = d2Jc[6 * 1 + 4] = -4 * s13 / 3;
		d2Jc[6 * 3 + 2] = d2Jc[6 * 2 + 3] = -4 * s12 / 3;
		d2Jc[6 * 3 + 0] = d2Jc[6 * 0 + 3] = +2 * s12 / 3;
		d2Jc[6 * 4 + 0] = d2Jc[6 * 0 + 4] = +2 * s13 / 3;
		d2Jc[6 * 3 + 1] = d2Jc[6 * 1 + 3] = +2 * s12 / 3;
		d2Jc[6 * 5 + 1] = d2Jc[6 * 1 + 5] = +2 * s23 / 3;
		d2Jc[6 * 4 + 2] = d2Jc[6 * 2 + 4] = +2 * s13 / 3;
		d2Jc[6 * 5 + 2] = d2Jc[6 * 2 + 5] = +2 * s23 / 3;
		d2Jc[6 * 1 + 0] = d2Jc[6 * 0 + 1] = 2 * (2 * s33 - s11 - s22) / 9;
		d2Jc[6 * 2 + 0] = d2Jc[6 * 0 + 2] = 2 * (2 * s22 - s11 - s33) / 9;
		d2Jc[6 * 2 + 1] = d2Jc[6 * 1 + 2] = 2 * (2 * s11 - s22 - s33) / 9;
		//return
		return stress_matrix(d2J3, d2Jc, t);
	}

	double von_mises_stress(const double* s, uint32_t t)
	{
		return sqrt(3 * deviatoric_invariant_J2(s, t));
	}
	double* von_mises_gradient(double* dse, const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//von misses stress
		const double se = von_mises_stress(s, t);
		//von misses gradient
		double dsc[6];
		dsc[3] = se != 0 ? 3 * s12 / se : 1;
		dsc[4] = se != 0 ? 3 * s13 / se : 1;
		dsc[5] = se != 0 ? 3 * s23 / se : 1;
		dsc[0] = se != 0 ? (s11 - s22 / 2 - s33 / 2) / se : 1;
		dsc[1] = se != 0 ? (s22 - s11 / 2 - s33 / 2) / se : 1;
		dsc[2] = se != 0 ? (s33 - s11 / 2 - s22 / 2) / se : 1;
		//return
		return stress_vector(dse, dsc, t);
	}
	double* von_mises_hessian(double* d2se, const double* s, uint32_t t)
	{
		//stress
		const double s11 = stress_value(s, t, (uint32_t) stress::s11);
		const double s22 = stress_value(s, t, (uint32_t) stress::s22);
		const double s33 = stress_value(s, t, (uint32_t) stress::s33);
		const double s12 = stress_value(s, t, (uint32_t) stress::s12);
		const double s13 = stress_value(s, t, (uint32_t) stress::s13);
		const double s23 = stress_value(s, t, (uint32_t) stress::s23);
		//von misses stress
		const double se = von_mises_stress(s, t);
		//von mises hessian
		double d2sc[36];
		d2sc[6 * 4 + 3] = d2sc[6 * 3 + 4] = -12 * s12 * s13;
		d2sc[6 * 5 + 3] = d2sc[6 * 3 + 5] = -12 * s12 * s23;
		d2sc[6 * 5 + 4] = d2sc[6 * 4 + 5] = -12 * s13 * s23;
		d2sc[6 * 3 + 0] = d2sc[6 * 0 + 3] = 2 * s12 * (s22 + s33 - 2 * s11);
		d2sc[6 * 4 + 0] = d2sc[6 * 0 + 4] = 2 * s13 * (s22 + s33 - 2 * s11);
		d2sc[6 * 5 + 0] = d2sc[6 * 0 + 5] = 2 * s23 * (s22 + s33 - 2 * s11);
		d2sc[6 * 3 + 1] = d2sc[6 * 1 + 3] = 2 * s12 * (s11 + s33 - 2 * s22);
		d2sc[6 * 4 + 1] = d2sc[6 * 1 + 4] = 2 * s13 * (s11 + s33 - 2 * s22);
		d2sc[6 * 5 + 1] = d2sc[6 * 1 + 5] = 2 * s23 * (s11 + s33 - 2 * s22);
		d2sc[6 * 3 + 2] = d2sc[6 * 2 + 3] = 2 * s12 * (s11 + s22 - 2 * s33);
		d2sc[6 * 4 + 2] = d2sc[6 * 2 + 4] = 2 * s13 * (s11 + s22 - 2 * s33);
		d2sc[6 * 5 + 2] = d2sc[6 * 2 + 5] = 2 * s23 * (s11 + s22 - 2 * s33);
		d2sc[6 * 0 + 0] = pow(s22 - s33, 2) + 4 * (s12 * s12 + s13 * s13 + s23 * s23);
		d2sc[6 * 1 + 1] = pow(s11 - s33, 2) + 4 * (s12 * s12 + s13 * s13 + s23 * s23);
		d2sc[6 * 2 + 2] = pow(s11 - s22, 2) + 4 * (s12 * s12 + s13 * s13 + s23 * s23);
		d2sc[6 * 1 + 0] = d2sc[6 * 0 + 1] = (s11 - s33) * (s33 - s22) - 2 * (s12 * s12 + s13 * s13 + s23 * s23);
		d2sc[6 * 2 + 0] = d2sc[6 * 0 + 2] = (s11 - s22) * (s22 - s33) - 2 * (s12 * s12 + s13 * s13 + s23 * s23);
		d2sc[6 * 2 + 1] = d2sc[6 * 1 + 2] = (s11 - s22) * (s33 - s11) - 2 * (s12 * s12 + s13 * s13 + s23 * s23);
		d2sc[6 * 3 + 3] = 4 * (s11 * s11 + s22 * s22 + s33 * s33 - s22 * s33 - s11 * (s22 + s33) + 3 * (s13 * s13 + s23 * s23));
		d2sc[6 * 4 + 4] = 4 * (s11 * s11 + s22 * s22 + s33 * s33 - s22 * s33 - s11 * (s22 + s33) + 3 * (s12 * s12 + s23 * s23));
		d2sc[6 * 5 + 5] = 4 * (s11 * s11 + s22 * s22 + s33 * s33 - s22 * s33 - s11 * (s22 + s33) + 3 * (s12 * s12 + s13 * s13));
		//return
		return stress_matrix(d2se, d2sc, t, se == 0 ? 1 : 4 * se * se * se / 3);
	}

	double lode_angle(const double* s, uint32_t t)
	{
		const double J2 = deviatoric_invariant_J2(s, t);
		const double J3 = deviatoric_invariant_J3(s, t);
		return fabs(J2) > 1e-5 ? acos(math::bound(3 * sqrt(3) / 2 * J3 / sqrt(J2 * J2 * J2))) / 3 : 0;
	}
}