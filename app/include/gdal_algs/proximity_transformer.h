#pragma once

#include "defs.h"
#include <gdal_priv.h>
#include "gdal_utils.h"
#include "gdal_alg.h"
#include "io_raster.h"


namespace LxGeo
{
	namespace lxProximityAlign
	{

		void transform_to_proximity(const std::string& in_raster_path, const std::string& out_raster_path);

	}
}