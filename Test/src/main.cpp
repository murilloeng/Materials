//std
#include <cstdlib>

//test
#include "Materials/Test/inc/test.hpp"
#include "Materials/Materials/inc/Mechanic/Stress.hpp"
#include "Materials/Materials/inc/Mechanic/Concrete.hpp"

int main(void)
{
	//data
	const bool steel = true;
	const bool stress = true;
	const bool concrete = true;
	const uint32_t st = uint32_t(materials::Stress::type::s11);
	const materials::Concrete::Model model = materials::Concrete::Model::bresler_pister;
	//test
	if(stress)
	{
		tests::stress::test_hessian_I2(st);
		tests::stress::test_hessian_I3(st);
		tests::stress::test_gradient_I1(st);
		tests::stress::test_gradient_I2(st);
		tests::stress::test_gradient_I3(st);
		tests::stress::test_von_mises_hessian(st);
		tests::stress::test_von_mises_gradient(st);
		tests::stress::test_hydrostatic_gradient(st);
		tests::stress::test_deviatoric_hessian_J2(st);
		tests::stress::test_deviatoric_hessian_J3(st);
		tests::stress::test_deviatoric_gradient_J2(st);
		tests::stress::test_deviatoric_gradient_J3(st);
	}
	if(steel)
	{
		tests::steel::test_uniaxial();
		tests::steel::test_yield_hessian(st);
		tests::steel::test_yield_gradient(st);
		tests::steel::test_plastic_gradient();
	}
	if(concrete)
	{
		tests::concrete::test_biaxial(model);
		tests::concrete::test_tension(model);
		tests::concrete::test_compression(model);
		tests::concrete::test_plastic_gradient();
		tests::concrete::test_yield_hessian(st, model);
		tests::concrete::test_yield_gradient(st, model);
	}
	//return
	return EXIT_SUCCESS;
}