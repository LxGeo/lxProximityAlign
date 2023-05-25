#pragma once
#include "defs.h"
#include "io_raster.h"
#include "geometries_with_attributes/geometries_with_attributes.h"
#include "parameters.h"
#include "io_shapefile.h"
#include "polygon_decompoe_utils.h"
#include "proximity_triplet.h"
#include "affine_geometry/affine_transformer.h"
#include "graph_weights/spatial_weights.h"
#include "relationships/composition_struct.h"


namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> alignmentNeighbour(std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster, std::vector<Geometries_with_attributes<Boost_Polygon_2>>& input_polygons);

	}
}