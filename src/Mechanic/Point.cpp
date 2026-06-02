//std
#include <cstring>

//math
#include "Math/inc/Miscellaneous/bits.hpp"

//Materials
#include "Materials/inc/Mechanic/Point.hpp"
#include "Materials/inc/Mechanic/States.hpp"

namespace materials
{
	//constructors
	Point::Point(void) :
		m_state_old(state::elastic), m_state_new(state::elastic), 
		m_stress_types(0), m_hardening_old(0), m_hardening_new(0), 
		m_plastic_strain_old(nullptr), m_plastic_strain_new(nullptr)
	{
		return;
	}

	//destructor
	Point::~Point(void)
	{
		delete[] m_plastic_strain_old;
		delete[] m_plastic_strain_new;
	}

	//analysis
	void Point::reset(void)
	{
		m_hardening_old = 0;
		m_hardening_new = 0;
		m_state_old = state::elastic;
		m_state_new = state::elastic;
		const uint32_t ns = math::bit_count(m_stress_types);
		memset(m_plastic_strain_old, 0, ns * sizeof(double));
		memset(m_plastic_strain_new, 0, ns * sizeof(double));
	}
	void Point::prepare(uint32_t stress_types)
	{
		//sizes
		m_stress_types = stress_types;
		const uint32_t ns = math::bit_count(m_stress_types);
		//state
		m_hardening_old = 0;
		m_hardening_new = 0;
		m_state_old = state::elastic;
		m_state_new = state::elastic;
		delete[] m_plastic_strain_old;
		delete[] m_plastic_strain_new;
		m_plastic_strain_old = new double[ns];
		m_plastic_strain_new = new double[ns];
		memset(m_plastic_strain_old, 0, ns * sizeof(double));
		memset(m_plastic_strain_new, 0, ns * sizeof(double));
	}

	void Point::update(void)
	{
		m_state_old = m_state_new;
		m_hardening_old = m_hardening_new;
		memcpy(m_plastic_strain_old, m_plastic_strain_new, math::bit_count(m_stress_types) * sizeof (double));
	}
	void Point::restore(void)
	{
		m_state_new = m_state_old;
		m_hardening_new = m_hardening_old;
		memcpy(m_plastic_strain_new, m_plastic_strain_old, math::bit_count(m_stress_types) * sizeof (double));
	}
}