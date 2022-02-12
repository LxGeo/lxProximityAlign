#pragma once
#include "defs.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{

		enum SupportPointsStrategy
		{
			vertex_only = 1 << 0,
			vertex_and_mid_point = 1 << 1,
			constant_walker = 1 << 2
		};

		struct SupportPoints {
			std::vector<size_t> polygon_indices;
			std::vector<Boost_Point_2> support_points;
		};

		SupportPoints decompose_polygons(std::vector<Boost_Polygon_2>& input_polygons,
			SupportPointsStrategy decompose_strategy = SupportPointsStrategy::vertex_and_mid_point);

	}
}