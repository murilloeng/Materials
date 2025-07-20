//materials
#include "Materials/Materials/inc/Mechanic/Steel.hpp"
#include "Materials/Materials/inc/Mechanic/Stress.hpp"

//test
#include "Materials/Test/inc/test.hpp"

void tests::steel::test_uniaxial(void)
{
	//data
	double e, s, C;
	materials::Steel steel;
	materials::Point point;
	const uint32_t np = 200;
	const double E = 2.00e+11;
	const double sy = 4.00e+08;
	const double Kp = 0.00e+00;
	materials::Mechanic::hardening(true);
	point.prepare(materials::Stress::uniaxial);
	//stream
	FILE* stream = popen("gnuplot.exe -p", "w");
	fprintf(stream, "set grid; plot '-' with lines notitle\n");
	//test
	steel.yield_stress(sy);
	steel.elastic_modulus(E);
	steel.plastic_modulus(Kp);
	for(uint32_t i = 0; i < np; i++)
	{
		e = 2 * sy / E * i / np;
		steel.return_mapping(&s, &C, &e, point);
		fprintf(stream, "%.6e %.6e\n", e, s);
	}
	fprintf(stream, "e\n");
	pclose(stream);
}
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