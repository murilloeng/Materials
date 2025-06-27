#pragma once

//std
#include <cstdint>

namespace materials
{
	enum class strain_measure : uint32_t
	{
		linear = 1 << 0,
		quadratic = 1 << 1,
		logarithmic = 1 << 2,
		last
	};
}