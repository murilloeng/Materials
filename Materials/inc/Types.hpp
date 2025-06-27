#pragma once

//std
#include <cstdint>

namespace materials
{
	enum class type : uint32_t
	{
		heat			= 1 << 0,
		steel			= 1 << 1,
		concrete		= 1 << 2,
		hyperelastic	= 1 << 3,
		last
	};
}