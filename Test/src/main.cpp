//std
#include <cstdlib>

//test
#include "Materials/Test/inc/test.hpp"
#include "Materials/Materials/inc/Mechanic/stress.hpp"

int main(void)
{
	//test
	const uint32_t st = 
		uint32_t(materials::stress::s11) |
		uint32_t(materials::stress::s22) |
		uint32_t(materials::stress::s33) |
		uint32_t(materials::stress::s12) |
		uint32_t(materials::stress::s13) |
		uint32_t(materials::stress::s23);
	tests::steel::test_yield_hessian(st);
	//return
	return EXIT_SUCCESS;
}