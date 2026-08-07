//std
#include <cstdlib>

//test
#include "Materials/Test/inc/test.hpp"
#include "Materials/inc/Mechanic/Stress.hpp"
#include "Materials/inc/Mechanic/Concrete.hpp"

int main(void)
{
	//data
	const bool steel = false;
	const bool stress = false;
	const bool concrete = true;
	const bool uniaxial = false;
	const uint32_t st = 1 << uint32_t(materials::Stress::Type::s11);
	const materials::Concrete::Model model = materials::Concrete::Model::bresler_pister;
	//test
	if(stress)
	{
		tests::stress::hessian_I2(st);
		tests::stress::hessian_I3(st);
		tests::stress::gradient_I1(st);
		tests::stress::gradient_I2(st);
		tests::stress::gradient_I3(st);
		tests::stress::von_mises_hessian(st);
		tests::stress::von_mises_gradient(st);
		tests::stress::hydrostatic_gradient(st);
		tests::stress::deviatoric_hessian_J2(st);
		tests::stress::deviatoric_hessian_J3(st);
		tests::stress::deviatoric_gradient_J2(st);
		tests::stress::deviatoric_gradient_J3(st);
	}
	if(steel)
	{
		tests::steel::uniaxial();
		tests::steel::yield_hessian(st);
		tests::steel::yield_gradient(st);
		tests::steel::plastic_gradient();
	}
	if(concrete)
	{
		tests::concrete::biaxial(model);
		tests::concrete::tension(model);
		tests::concrete::compression(model);
		tests::concrete::plastic_gradient();
		tests::concrete::yield_hessian(st, model);
		tests::concrete::yield_gradient(st, model);
	}
	if(uniaxial)
	{
		tests::uniaxial::tension();
		tests::uniaxial::compression();
	}
	//return
	return EXIT_SUCCESS;
}