#pragma once
#include "defs.h"
#include "io_shapefile.h"
#include <gdal_alg.h>
#include "io_raster.h"
#include "gdal_algs/rasterizer.h"
#include "gdal_algs/proximity_transformer.h"
#include "parameters.h"
#include "geometry_lab.h"

namespace LxGeo
{
    namespace lxProximityAlign
    {

		enum class ProximityMapStrategy
		{
			vertex_only = 1 << 0,
			contours = 1 << 1,
			filled_polygon = 1 << 2,
			constant_walker = 1 << 3,
			skeleton =  1 << 4
		};

		void polygons2proximity(IO_DATA::PolygonsShapfileIO& input_shapefile, std::string& out_proximity_map_path, OGREnvelope* out_extents =NULL,
			double raster_px_size=0.5, double raster_py_size=0.5, ProximityMapStrategy startegy = ProximityMapStrategy::contours);

    }
}