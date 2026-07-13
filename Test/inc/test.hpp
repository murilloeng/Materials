#pragma once

//std
#include <cstdint>

//Materials
#include "Materials/inc/Mechanic/Concrete.hpp"

namespace tests
{
	namespace stress
	{
		void hessian_I2(uint32_t);
		void hessian_I3(uint32_t);

		void gradient_I1(uint32_t);
		void gradient_I2(uint32_t);
		void gradient_I3(uint32_t);

		void von_mises_hessian(uint32_t);
		void von_mises_gradient(uint32_t);

		void hydrostatic_gradient(uint32_t);

		void deviatoric_hessian_J2(uint32_t);
		void deviatoric_hessian_J3(uint32_t);

		void deviatoric_gradient_J2(uint32_t);
		void deviatoric_gradient_J3(uint32_t);
	}
	namespace steel
	{
		void uniaxial(void);
		void plastic_gradient(void);
		void yield_hessian(uint32_t);
		void yield_gradient(uint32_t);
	}
	namespace concrete
	{
		void plastic_gradient(void);
		void tension(materials::Concrete::Model);
		void biaxial(materials::Concrete::Model);
		void compression(materials::Concrete::Model);
		void yield_hessian(uint32_t, materials::Concrete::Model);
		void yield_gradient(uint32_t, materials::Concrete::Model);
	}
	namespace uniaxial
	{
		void tension(void);
		void compression(void);
	}
}