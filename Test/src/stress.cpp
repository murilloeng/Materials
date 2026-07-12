//std
#include <cmath>
#include <ctime>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Miscellaneous/bits.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

//Materials
#include "Materials/inc/Mechanic/Stress.hpp"

//test
#include "Materials/Test/inc/test.hpp"

void tests::stress::test_hessian_I2(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Matrix s(ns, 1), ddI2a(ns, ns), ddI2n(ns, ns), ddI2r(ns, ns);
	std::function<void(double*, const double*)> function = [st](double* dI2, const double* s){
		materials::Stress(s, st).gradient_I2(dI2);
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		materials::Stress(s.data(), st).hessian_I2(ddI2a.data());
		math::ndiff(function, ddI2n.data(), s.data(), ns, ns, 1e-7);
		ddI2r = ddI2a - ddI2n;
		if(ddI2r.norm() > 1e-5 * ddI2a.norm())
		{
			s.print("s");
			ddI2a.print("dI2a");
			ddI2n.print("dI2n");
			ddI2r.print("dI2r", 1e-5);
			printf("Hessian I2 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Hessian I2 test: %d ok!\n", i);
		}
	}
}
void tests::stress::test_hessian_I3(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Matrix s(ns, 1), ddI3a(ns, ns), ddI3n(ns, ns), ddI3r(ns, ns);
	std::function<void(double*, const double*)> function = [st](double* dI3, const double* s){
		materials::Stress(s, st).gradient_I3(dI3);
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		materials::Stress(s.data(), st).hessian_I3(ddI3a.data());
		math::ndiff(function, ddI3n.data(), s.data(), ns, ns, 1e-7);
		ddI3r = ddI3a - ddI3n;
		if(ddI3r.norm() > 1e-5 * ddI3a.norm())
		{
			s.print("s");
			ddI3a.print("dI3a");
			ddI3n.print("dI3n");
			ddI3r.print("dI3r", 1e-5);
			printf("Hessian I3 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Hessian I3 test: %d ok!\n", i);
		}
	}
}

void tests::stress::test_gradient_I1(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Vector s(ns), dI1a(ns), dI1n(ns), dI1r(ns);
	std::function<void(double*, const double*)> function = [st](double* I1, const double* s){
		I1[0] = materials::Stress(s, st).invariant_I1();
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		materials::Stress(s.data(), st).gradient_I1(dI1a.data());
		math::ndiff(function, dI1n.data(), s.data(), 1, ns, 1e-7);
		dI1r = dI1a - dI1n;
		if(dI1r.norm() > 1e-5 * dI1a.norm())
		{
			s.print("s");
			dI1a.print("dI1a");
			dI1n.print("dI1n");
			dI1r.print("dI1r", 1e-5);
			printf("Gradient I1 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Gradient I1 test: %d ok!\n", i);
		}
	}
}
void tests::stress::test_gradient_I2(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Vector s(ns), dI2a(ns), dI2n(ns), dI2r(ns);
	std::function<void(double*, const double*)> function = [st](double* I2, const double* s){
		I2[0] = materials::Stress(s, st).invariant_I2();
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		materials::Stress(s.data(), st).gradient_I2(dI2a.data());
		math::ndiff(function, dI2n.data(), s.data(), 1, ns, 1e-7);
		dI2r = dI2a - dI2n;
		if(dI2r.norm() > 1e-5 * dI2a.norm())
		{
			s.print("s");
			dI2a.print("dI2a");
			dI2n.print("dI2n");
			dI2r.print("dI2r", 1e-5);
			printf("Gradient I2 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Gradient I2 test: %d ok!\n", i);
		}
	}
}
void tests::stress::test_gradient_I3(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Vector s(ns), dI3a(ns), dI3n(ns), dI3r(ns);
	std::function<void(double*, const double*)> function = [st](double* I3, const double* s){
		I3[0] = materials::Stress(s, st).invariant_I3();
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		materials::Stress(s.data(), st).gradient_I3(dI3a.data());
		math::ndiff(function, dI3n.data(), s.data(), 1, ns, 1e-7);
		dI3r = dI3a - dI3n;
		if(dI3r.norm() > 1e-5 * dI3a.norm())
		{
			s.print("s");
			dI3a.print("dI3a");
			dI3n.print("dI3n");
			dI3r.print("dI3r", 1e-5);
			printf("Gradient I3 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Gradient I3 test: %d ok!\n", i);
		}
	}
}

void tests::stress::test_von_mises_hessian(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Matrix s(ns, 1), ddsea(ns, ns), ddsen(ns, ns), ddser(ns, ns);
	std::function<void(double*, const double*)> function = [st](double* dse, const double* s){
		materials::Stress(s, st).von_mises_gradient(dse);
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		math::ndiff(function, ddsen.data(), s.data(), ns, ns, 1e-7);
		materials::Stress(s.data(), st).von_mises_hessian(ddsea.data());
		ddser = ddsea - ddsen;
		if(ddser.norm() > 1e-5 * (ddsea.norm() ? ddsea.norm() : 1))
		{
			s.print("s");
			ddsea.print("ddsea");
			ddsen.print("ddsen");
			ddser.print("ddser", 1e-5);
			printf("von Mises hessian numerical test failed!\n");
			break;
		}
		else
		{
			printf("von Mises hessian test: %d ok!\n", i);
		}
	}
}
void tests::stress::test_von_mises_gradient(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Vector s(ns), dsea(ns), dsen(ns), dser(ns);
	std::function<void(double*, const double*)> function = [st](double* se, const double* s){
		se[0] = materials::Stress(s, st).von_mises_stress();
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		math::ndiff(function, dsen.data(), s.data(), 1, ns, 1e-7);
		materials::Stress(s.data(), st).von_mises_gradient(dsea.data());
		dser = dsea - dsen;
		if(dser.norm() > 1e-5 * dsea.norm())
		{
			s.print("s");
			dsea.print("dsea");
			dsen.print("dsen");
			dser.print("dser", 1e-5);
			printf("von Mises gradient numerical test failed!\n");
			break;
		}
		else
		{
			printf("von Mises gradient test: %d ok!\n", i);
		}
	}
}

void tests::stress::test_hydrostatic_gradient(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Vector s(ns), dsma(ns), dsmn(ns), dsmr(ns);
	std::function<void(double*, const double*)> function = [st](double* sm, const double* s){
		sm[0] = materials::Stress(s, st).hydrostatic_stress();
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		math::ndiff(function, dsmn.data(), s.data(), 1, ns, 1e-7);
		materials::Stress(s.data(), st).hydrostatic_gradient(dsma.data());
		dsmr = dsma - dsmn;
		if(dsmr.norm() > 1e-5 * dsma.norm())
		{
			s.print("s");
			dsma.print("dsma");
			dsmn.print("dsmn");
			dsmr.print("dsmr", 1e-5);
			printf("Hydrostatic stress gradient numerical test failed!\n");
			break;
		}
		else
		{
			printf("Hydrostatic stress gradient test: %d ok!\n", i);
		}
	}
}

void tests::stress::test_deviatoric_hessian_J2(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Matrix s(ns, 1), ddJ2a(ns, ns), ddJ2n(ns, ns), ddJ2r(ns, ns);
	std::function<void(double*, const double*)> function = [st](double* dJ2, const double* s){
		materials::Stress(s, st).deviatoric_gradient_J2(dJ2);
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		math::ndiff(function, ddJ2n.data(), s.data(), ns, ns, 1e-7);
		materials::Stress(s.data(), st).deviatoric_hessian_J2(ddJ2a.data());
		ddJ2r = ddJ2a - ddJ2n;
		if(ddJ2r.norm() > 1e-5 * ddJ2a.norm())
		{
			s.print("s");
			ddJ2a.print("ddJ2a");
			ddJ2n.print("ddJ2n");
			ddJ2r.print("ddJ2r", 1e-5);
			printf("Deviatoric stress hessian J2 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Deviatoric stress hessian J2 test: %d ok!\n", i);
		}
	}
}
void tests::stress::test_deviatoric_hessian_J3(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Matrix s(ns, 1), ddJ3a(ns, ns), ddJ3n(ns, ns), ddJ3r(ns, ns);
	std::function<void(double*, const double*)> function = [st](double* dJ3, const double* s){
		materials::Stress(s, st).deviatoric_gradient_J3(dJ3);
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		math::ndiff(function, ddJ3n.data(), s.data(), ns, ns, 1e-7);
		materials::Stress(s.data(), st).deviatoric_hessian_J3(ddJ3a.data());
		ddJ3r = ddJ3a - ddJ3n;
		if(ddJ3r.norm() > 1e-5 * ddJ3a.norm())
		{
			s.print("s");
			ddJ3a.print("ddJ3a");
			ddJ3n.print("ddJ3n");
			ddJ3r.print("ddJ3r", 1e-5);
			printf("Deviatoric stress hessian J3 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Deviatoric stress hessian J3 test: %d ok!\n", i);
		}
	}
}

void tests::stress::test_deviatoric_gradient_J2(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Vector s(ns), dJ2a(ns), dJ2n(ns), dJ2r(ns);
	std::function<void(double*, const double*)> function = [st](double* J2, const double* s){
		J2[0] = materials::Stress(s, st).deviatoric_invariant_J2();
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		math::ndiff(function, dJ2n.data(), s.data(), 1, ns, 1e-7);
		materials::Stress(s.data(), st).deviatoric_gradient_J2(dJ2a.data());
		dJ2r = dJ2a - dJ2n;
		if(dJ2r.norm() > 1e-5 * dJ2a.norm())
		{
			s.print("s");
			dJ2a.print("dJ2a");
			dJ2n.print("dJ2n");
			dJ2r.print("dJ2r", 1e-5);
			printf("Deviatoric stress gradient J2 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Deviatoric stress gradient J2 test: %d ok!\n", i);
		}
	}
}
void tests::stress::test_deviatoric_gradient_J3(uint32_t st)
{
	//data
	const uint32_t nt = 100000;
	const uint32_t ns = math::bit_count(st);
	math::Vector s(ns), dJ3a(ns), dJ3n(ns), dJ3r(ns);
	std::function<void(double*, const double*)> function = [st](double* J3, const double* s){
		J3[0] = materials::Stress(s, st).deviatoric_invariant_J3();
	};
	//test
	srand(uint32_t(time(nullptr)));
	for(uint32_t i = 0; i < nt; i++)
	{
		s.randu();
		math::ndiff(function, dJ3n.data(), s.data(), 1, ns, 1e-7);
		materials::Stress(s.data(), st).deviatoric_gradient_J3(dJ3a.data());
		dJ3r = dJ3a - dJ3n;
		if(dJ3r.norm() > 1e-5 * dJ3a.norm())
		{
			s.print("s");
			dJ3a.print("dJ3a");
			dJ3n.print("dJ3n");
			dJ3r.print("dJ3r", 1e-5 * dJ3a.norm());
			printf("Deviatoric stress gradient J3 numerical test failed!\n");
			break;
		}
		else
		{
			printf("Deviatoric stress gradient J3 test: %d ok!\n", i);
		}
	}
}