//materials
#include "Materials/Materials/inc/Mechanic/Steel.hpp"

//test
#include "Materials/Test/inc/test.hpp"

void tests::steel::test_plastic_gradient(void)
{
	materials::Steel steel;
	steel.plastic_modulus(1e5);
	steel.test_plastic_gradient();
}
void tests::steel::test_yield_hessian(uint32_t st)
{
	materials::Steel steel;
	steel.test_yield_hessian(st);
}
void tests::steel::test_yield_gradient(uint32_t st)
{
	materials::Steel steel;
	steel.test_yield_gradient(st);
}