#pragma once
#include "defs.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{

		/**
		*  A DhmTransformer class to manage running required steps to generate final transformed raster.
		*/
		class OptimProximityAligner
		{

		public:
			OptimProximityAligner() {};

			~OptimProximityAligner() {};

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

			double couple_rotation_angle = DBL_MAX;
			double couple_v_displacement = 0;
			std::pair<double, double> xy_cst;

		};
	}
}