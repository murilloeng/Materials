#pragma once

//std
#include <cstdint>

namespace materials
{
	enum class state : uint32_t
	{
		broken,
		elastic,
		plastic
	};
}