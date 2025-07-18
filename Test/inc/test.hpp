#pragma once

//std
#include <cstdint>

namespace tests
{
	namespace steel
	{
		void test_plastic_gradient(void);
		void test_yield_hessian(uint32_t);
		void test_yield_gradient(uint32_t);
	}
}