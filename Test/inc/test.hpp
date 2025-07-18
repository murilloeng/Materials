#pragma once

//std
#include <cstdint>

//materials
#include "Materials/Materials/inc/Mechanic/Concrete.hpp"

namespace tests
{
	namespace stress
	{
		void test_hessian_I2(uint32_t);
		void test_hessian_I3(uint32_t);

		void test_gradient_I1(uint32_t);
		void test_gradient_I2(uint32_t);
		void test_gradient_I3(uint32_t);

		void test_von_mises_hessian(uint32_t);
		void test_von_mises_gradient(uint32_t);

		void test_hydrostatic_gradient(uint32_t);

		void test_deviatoric_hessian_J2(uint32_t);
		void test_deviatoric_hessian_J3(uint32_t);

		void test_deviatoric_gradient_J2(uint32_t);
		void test_deviatoric_gradient_J3(uint32_t);
	}
	namespace steel
	{
		void test_uniaxial(void);
		void test_plastic_gradient(void);
		void test_yield_hessian(uint32_t);
		void test_yield_gradient(uint32_t);
	}
	namespace concrete
	{
		void test_plastic_gradient(void);
		void test_yield_hessian(uint32_t, materials::Concrete::Model);
		void test_yield_gradient(uint32_t, materials::Concrete::Model);
	}
}