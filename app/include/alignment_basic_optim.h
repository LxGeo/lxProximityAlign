#pragma once
#include "defs.h"
#include "io_raster.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> alignmentIteration(std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster, std::vector<Boost_Polygon_2>& input_polygons);

	}
}