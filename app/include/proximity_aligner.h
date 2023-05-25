#pragma once
#include "defs.h"
#include "io_raster.h"
#include "io_shapefile.h"
#include "defs.h"
#include "parameters.h"
#include "graph_weights/spatial_weights.h"
#include "geometries_with_attributes/linestring_with_attributes.h"
//#include "alignment_basic_optim.h"
//#include "alignment_iterative_support_optim.h"
#include "alignment_neighbour_optim.h"
//#include "alignment_neighbour_proximity_aware.h"

#include "gdal_algs/polygons_to_proximity_map.h"
#include "lightweight/geovector.h"


namespace LxGeo
{
	namespace lxProximityAlign
	{

		/**
		*  A DhmTransformer class to manage running required steps to generate final transformed raster.
		*/
		class ProximityAligner
		{

		public:
			ProximityAligner() {};

			~ProximityAligner() {};

			/**
			*  A method used to run all steps of transformation.
			*/
			virtual void run();

			/**
			*  A method to check requirements before running transformation steps.
			* Example: -Checking input_raster exsitance, checking output_path overwrite, check algorithm parameters ...
			* @return an bool indicating if can run algorithm
			*/
			bool pre_check();

		};
	}
}

