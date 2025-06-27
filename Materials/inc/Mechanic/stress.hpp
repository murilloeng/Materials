#pragma once

//std
#include <cstdint>

namespace materials
{
	enum class stress : uint32_t
	{
		s11 = 1 << 0,
		s22 = 1 << 1,
		s33 = 1 << 2,
		s12 = 1 << 3,
		s13 = 1 << 4,
		s23 = 1 << 5,
		last
	};

	double stress_value(const double*, uint32_t, uint32_t);
	double* stress_vector(double*, const double*, uint32_t, double = 1);
	double* stress_matrix(double*, const double*, uint32_t, double = 1);

	double stress_invariant_I1(const double*, uint32_t);
	double stress_invariant_I2(const double*, uint32_t);
	double stress_invariant_I3(const double*, uint32_t);
	double stress_principal_s1(const double*, uint32_t);
	double stress_principal_s2(const double*, uint32_t);
	double stress_principal_s3(const double*, uint32_t);

	double* stress_gradient_I1(double*, const double*, uint32_t);
	double* stress_gradient_I2(double*, const double*, uint32_t);
	double* stress_gradient_I3(double*, const double*, uint32_t);

	double* stress_hessian_I2(double*, const double*, uint32_t);
	double* stress_hessian_I3(double*, const double*, uint32_t);

	double hydrostatic_stress(const double*, uint32_t);
	double* hydrostatic_gradient(double*, const double*, uint32_t);

	double deviatoric_invariant_J2(const double*, uint32_t);
	double deviatoric_invariant_J3(const double*, uint32_t);
	double deviatoric_principal_s1(const double*, uint32_t);
	double deviatoric_principal_s2(const double*, uint32_t);
	double deviatoric_principal_s3(const double*, uint32_t);

	double* deviatoric_gradient_J2(double*, const double*, uint32_t);
	double* deviatoric_gradient_J3(double*, const double*, uint32_t);

	double* deviatoric_hessian_J2(double*, const double*, uint32_t);
	double* deviatoric_hessian_J3(double*, const double*, uint32_t);

	double von_mises_stress(const double*, uint32_t);
	double* von_mises_gradient(double*, const double*, uint32_t);
	double* von_mises_hessian(double*, const double*, uint32_t);

	double lode_angle(const double*, uint32_t);
}