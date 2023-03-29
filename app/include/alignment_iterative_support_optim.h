#pragma once
#include "defs.h"
#include "io_raster.h"
#include "geometries_with_attributes/geometries_with_attributes.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> iterative_support_alignment(std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster, std::vector<Geometries_with_attributes<Boost_Polygon_2>>& input_polygons);

	}
}