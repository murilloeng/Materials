#pragma once

//std
#include <cstdint>

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			enum class hyperelastic_type : uint32_t
			{
				gent,
				hencky,
				neo_hookean,
				saint_venant,
				last
			};
		}
	}
}