//std
#include <cmath>

//Math
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Validation/Validator.hpp"

//Materials
#include "Materials/inc/Mechanic/Point.hpp"
#include "Materials/inc/Mechanic/Stress.hpp"
#include "Materials/inc/Mechanic/Uniaxial.hpp"

//Test
#include "Materials/Test/inc/test.hpp"

//data
static double s, C;
static const uint32_t np = 100;
static const double E = 2.00e+11;
static const double K = 1.00e+11;
static const double sy = 4.00e+08;
static const double sr = 2.00e+08;

static double data[2 * (np + 1)];

static double function(double e)
{
	const double Et = E * K / (E + K);
	const double ey = (sy - sr * math::sign(e)) / E;
	return fabs(E * e + sr) < sy ? E * e + sr : sy * math::sign(e) + Et * (e - ey * math::sign(e));
}

void tests::uniaxial::tension(void)
{
	//data
	materials::Point point;
	materials::Uniaxial material;
	math::validation::Validator validator;
	const uint32_t st = uint32_t(materials::Stress::type::s11);
	//return mapping
	point.prepare(st);
	material.inelastic(true);
	material.yield_stress(sy);
	material.elastic_modulus(E);
	material.plastic_modulus(K);
	for(uint32_t i = 0; i <= np; i++)
	{
		const double e = +2 * sy / E * i / np;
		material.return_mapping(s, C, e, sr, point);
		data[0 * (np + 1) + i] = e, data[1 * (np + 1) + i] = s;
		printf("Strain: %+.2e Stress: %+.2e Stiffness: %+.2e State: %d\n", e, s, C, uint32_t(point.m_state_new));
		point.update();
	}
	//validatior
	validator.create_item();
	validator.item(0)->function(function);
	validator.item(0)->load_numeric(data, data + np + 1, np + 1);
	//validation
	validator.validate();
}
void tests::uniaxial::compression(void)
{
	//data
	materials::Point point;
	materials::Uniaxial material;
	math::validation::Validator validator;
	const uint32_t st = uint32_t(materials::Stress::type::s11);
	//return mapping
	point.prepare(st);
	material.inelastic(true);
	material.yield_stress(sy);
	material.elastic_modulus(E);
	material.plastic_modulus(K);
	for(uint32_t i = 0; i <= np; i++)
	{
		const double e = -2 * sy / E * i / np;
		material.return_mapping(s, C, e, sr, point);
		data[0 * (np + 1) + i] = e, data[1 * (np + 1) + i] = s;
		printf("Strain: %+.2e Stress: %+.2e Stiffness: %+.2e State: %d\n", e, s, C, uint32_t(point.m_state_new));
		point.update();
	}
	//validatior
	validator.create_item();
	validator.item(0)->function(function);
	validator.item(0)->load_numeric(data, data + np + 1, np + 1);
	//validation
	validator.validate();
}