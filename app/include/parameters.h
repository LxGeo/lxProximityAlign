#pragma once
#include "defs.h"
#include "cli/base_parameters.h"
#include "CLI/CLI.hpp"
#include <boost/filesystem.hpp>


namespace LxGeo
{
	namespace lxProximityAlign
	{
		class Parameters
		{
		public:
			Parameters(int argc, char* argv[]);

			~Parameters();

			bool initialized();

		protected:
			void init();

		public:
			bool printed_help;

			std::string input_shapefile_to_align;
			std::string input_ref_shapefile;

			std::string output_basename;
			std::string output_shapefile;

			std::string temp_dir;

			bool overwrite_output;
			bool keep_geometries;

			double max_disparity;

			std::string r_imd_path;
			std::string i_imd_path;

		};

		extern Parameters* params;
	}
}
