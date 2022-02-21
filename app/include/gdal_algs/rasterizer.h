#pragma once
#include "defs.h"


namespace LxGeo
{
	namespace lxProximityAlign
	{

		void rasterize_shapefile(const std::string& out_raster_path, const std::string& input_shapefile_path, OGREnvelope* raster_extents = NULL, double raster_px_size = 0.5, double raster_py_size = 0.5);

	}
}