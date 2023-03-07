#pragma once
#include "defs.h"
#include "io_raster.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> nm_proximity_align(std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster, std::vector<Boost_Polygon_2>& input_polygons);

		std::vector<Boost_Polygon_2> nm_proximity_align_1d(
			std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster,
			std::vector<Boost_Polygon_2>& input_polygons,
			std::pair<double, double>& r2r_constants
		);

	}
}