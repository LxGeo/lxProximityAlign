#pragma once
#include "defs.h"

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

