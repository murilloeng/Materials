#pragma once

//std
#include <cstdint>

namespace materials
{
	enum class Type : uint32_t
	{
		heat			= 1 << 0,
		steel			= 1 << 1,
		concrete		= 1 << 2,
		last
	};
}