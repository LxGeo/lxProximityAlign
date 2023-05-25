#pragma once
#include "defs.h"
#include "io_raster.h"
#include "parameters.h"
#include "io_shapefile.h"
#include "alignment_neighbour_optim.h"
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

		std::vector<Boost_Polygon_2> alignmentNeighbourProxAware(std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster, std::vector<Boost_Polygon_2>& input_polygons);

	}
}