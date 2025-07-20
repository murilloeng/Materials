//materials
#include "Materials/Materials/inc/Mechanic/Stress.hpp"
#include "Materials/Materials/inc/Mechanic/Concrete.hpp"

//test
#include "Materials/Test/inc/test.hpp"

void tests::concrete::test_tension(materials::Concrete::Model model)
{
	//data
	double e, s, C;
	materials::Point point;
	const uint32_t np = 200;
	materials::Concrete concrete;
	materials::Mechanic::hardening(true);
	point.prepare(materials::Stress::uniaxial);
	const double E = concrete.elastic_modulus();
	const double st = concrete.yield_stress_tension();
	//stream
	FILE* stream = popen("gnuplot.exe -p", "w");
	fprintf(stream, "set grid; plot '-' with lines notitle\n");
	//test
	concrete.model(model);
	for(uint32_t i = 0; i < np; i++)
	{
		e = 2 * st / E * i / np;
		concrete.return_mapping(&s, &C, &e, point);
		fprintf(stream, "%.6e %.6e\n", e, s);
	}
	fprintf(stream, "e\n");
	pclose(stream);
}
void tests::concrete::test_biaxial(materials::Concrete::Model model)
{
	return;
}
void tests::concrete::test_compression(materials::Concrete::Model model)
{
	//data
	double e, s, C;
	materials::Point point;
	const uint32_t np = 200;
	materials::Concrete concrete;
	materials::Mechanic::hardening(true);
	point.prepare(materials::Stress::uniaxial);
	const double E = concrete.elastic_modulus();
	const double sc = concrete.yield_stress_compression();
	//stream
	FILE* stream = popen("gnuplot.exe -p", "w");
	fprintf(stream, "set grid; plot '-' with lines notitle\n");
	//test
	concrete.model(model);
	for(uint32_t i = 0; i < np; i++)
	{
		e = -2 * sc / E * i / np;
		concrete.return_mapping(&s, &C, &e, point);
		fprintf(stream, "%.6e %.6e\n", e, s);
	}
	fprintf(stream, "e\n");
	pclose(stream);
}
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