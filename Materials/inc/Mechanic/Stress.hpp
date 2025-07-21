#pragma once

//std
#include <cstdint>

namespace materials
{
	class Stress
	{
	public:
		//constructor
		Stress(const double*, uint32_t);

		//destructor
		~Stress(void);

		//enums
		enum class type : uint32_t
		{
			s11 = 1 << 0, s22 = 1 << 1, s33 = 1 << 2,
			s12 = 1 << 3, s13 = 1 << 4, s23 = 1 << 5, last
		};

		//operators
		double operator()(type) const;
		double operator[](type) const;

		//invariants
		double invariant_I1(void) const;
		double invariant_I2(void) const;
		double invariant_I3(void) const;

		//principal
		double principal_s1(void) const;
		double principal_s2(void) const;
		double principal_s3(void) const;

		//gradients
		double* gradient_I1(double*) const;
		double* gradient_I2(double*) const;
		double* gradient_I3(double*) const;

		//hessians
		double* hessian_I2(double*) const;
		double* hessian_I3(double*) const;

		//hydrostatic
		double hydrostatic_stress(void) const;
		double* hydrostatic_gradient(double*) const;

		//deviatoric
		double deviatoric_invariant_J2(void) const;
		double deviatoric_invariant_J3(void) const;

		double deviatoric_principal_s1(void) const;
		double deviatoric_principal_s2(void) const;
		double deviatoric_principal_s3(void) const;

		double* deviatoric_gradient_J2(double*) const;
		double* deviatoric_gradient_J3(double*) const;

		double* deviatoric_hessian_J2(double*) const;
		double* deviatoric_hessian_J3(double*) const;

		//von mises
		double von_mises_stress(void) const;
		double* von_mises_hessian(double*) const;
		double* von_mises_gradient(double*) const;

		//angle
		double lode_angle(void) const;

		//filters
		double* filter_vector(double*, const double*, double = 1) const;
		double* filter_matrix(double*, const double*, double = 1) const;

	private:
		//data
		uint32_t m_types;
		const double* m_data;
	};
}