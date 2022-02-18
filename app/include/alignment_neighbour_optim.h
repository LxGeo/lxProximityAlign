#pragma once
#include "defs.h"
#include "io_raster.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> alignmentNeighbour(std::map<std::string, RasterIO>& rasters_map, std::vector<Boost_Polygon_2>& input_polygons);

	}
}