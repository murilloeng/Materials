//Materials
#include "Materials/inc/Mechanic/Stress.hpp"
#include "Materials/inc/Mechanic/Concrete.hpp"

//test
#include "Materials/Test/inc/test.hpp"

void tests::concrete::tension(materials::Concrete::Model model)
{
	//data
	double e, s, C;
	materials::Point point;
	const uint32_t np = 200;
	materials::Concrete concrete;
	materials::Mechanic::hardening(true);
	const double E = concrete.elastic_modulus();
	const double st = concrete.yield_stress_tension();
	point.setup(uint32_t(materials::Stress::Type::s11));
	//stream
	FILE* stream = popen("gnuplot -p", "w");
	fprintf(stream, "set grid; plot '-' with lines linecolor rgb '#0000ff' notitle\n");
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
void tests::concrete::biaxial(materials::Concrete::Model model)
{
	//data
	materials::Point point;
	const uint32_t np = 200;
	double e[2], s[2], C[4];
	materials::Concrete concrete;
	materials::Mechanic::hardening(true);
	const double E = concrete.elastic_modulus();
	const double sb = concrete.yield_stress_biaxial();
	point.setup(uint32_t(materials::Stress::Type::s11) | uint32_t(materials::Stress::Type::s22));
	//stream
	FILE* stream = popen("gnuplot -p", "w");
	fprintf(stream, "set grid; plot '-' with lines linecolor rgb '#0000ff' notitle\n");
	//test
	concrete.model(model);
	for(uint32_t i = 0; i < np; i++)
	{
		e[0] = -2 * sb / E * i / np;
		e[1] = -2 * sb / E * i / np;
		concrete.return_mapping(s, C, e, point);
		fprintf(stream, "%.6e %.6e\n", e[0], s[0]);
	}
	fprintf(stream, "e\n");
	pclose(stream);
}
void tests::concrete::compression(materials::Concrete::Model model)
{
	//data
	double e, s, C;
	materials::Point point;
	const uint32_t np = 200;
	materials::Concrete concrete;
	materials::Mechanic::hardening(true);
	const double E = concrete.elastic_modulus();
	const double sc = concrete.yield_stress_compression();
	point.setup(uint32_t(materials::Stress::Type::s11));
	//stream
	FILE* stream = popen("gnuplot -p", "w");
	fprintf(stream, "set grid; plot '-' with lines linecolor rgb '#0000ff' notitle\n");
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
void tests::concrete::plastic_gradient(void)
{
	materials::Concrete concrete;
	concrete.test_plastic_gradient();
}
void tests::concrete::yield_hessian(uint32_t st, materials::Concrete::Model model)
{
	materials::Concrete concrete;
	concrete.model(model);
	concrete.test_yield_hessian(st);
}
void tests::concrete::yield_gradient(uint32_t st, materials::Concrete::Model model)
{
	materials::Concrete concrete;
	concrete.model(model);
	concrete.test_yield_gradient(st);
}