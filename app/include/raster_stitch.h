#pragma once
#include "defs.h"
#include "io_raster.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{
		enum class RasterPixelsStitcherStartegy {
			vertex_only = 1 << 0,
			vertex_and_mid_point = 1 << 1,
			constant_walker = 1 << 2,
			contours = 1 << 3
		};


		using namespace LxGeo::IO_DATA;

		class RasterPixelsStitcher {

		public:
			RasterPixelsStitcher() {};
			RasterPixelsStitcher(RasterIO& _ref_raster) { ref_raster = _ref_raster; };

			double readPolygonsPixels(std::list<Boost_Polygon_2>& resp_polygons, RasterPixelsStitcherStartegy strategy);

		public:
			RasterIO ref_raster;

		};
	}
}