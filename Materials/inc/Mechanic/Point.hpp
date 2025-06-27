#pragma once

//std
#include <cstdint>

namespace materials
{
	enum class state : uint32_t;
}

namespace materials
{
	class Point
	{
	public:
		//constructors
		Point(void);

		//destructor
		~Point(void);

		//analysis
		void reset(void);
		void prepare(uint32_t);

		void update(void);
		void restore(void);

		//data
		state m_state_old;
		state m_state_new;

		uint32_t m_stress_types;

		double m_hardening_old;
		double m_hardening_new;
		double* m_plastic_strain_old;
		double* m_plastic_strain_new;
	};
}