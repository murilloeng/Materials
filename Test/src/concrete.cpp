//materials
#include "Materials/Materials/inc/Mechanic/Concrete.hpp"

//test
#include "Materials/Test/inc/test.hpp"

void tests::concrete::test_plastic_gradient(void)
{
	materials::Concrete concrete;
	concrete.test_plastic_gradient();
}
void tests::concrete::test_yield_hessian(uint32_t st, materials::Concrete::Model model)
{
	materials::Concrete concrete;
	concrete.model(model);
	concrete.test_yield_hessian(st);
}
void tests::concrete::test_yield_gradient(uint32_t st, materials::Concrete::Model model)
{
	materials::Concrete concrete;
	concrete.model(model);
	concrete.test_yield_gradient(st);
}