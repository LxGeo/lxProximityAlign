#pragma once
#include "defs.h"
#include "io_raster.h"
#include "io_shapefile.h"
#include "parameters.h"
#include "lightweight/geovector.h"
#include "gdal_algs/polygons_to_proximity_map.h"
#include "nm_search.h"
#include "lightweight/geovector.h"
#include "affine_geometry/affine_transformer.h"
#include "stitching/vector_on_raster_stitcher.h"
#include "design_pattern/extended_iterators.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{

		/**
		*  A DhmTransformer class to manage running required steps to generate final transformed raster.
		*/
		class OptimProximityAlignerLinear
		{

		public:
			OptimProximityAlignerLinear() {};

			~OptimProximityAlignerLinear() {};

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

			double MAX_DISP = 100;

		};
	}
}