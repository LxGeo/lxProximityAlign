#pragma once
#include "defs.h"
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

			void help();

			void parse(int argc, char* argv[]);

		public:
			bool printed_help;

			std::string input_shapefile_to_align;
			std::string input_ref_shapefile;

			std::string output_basename;
			std::string output_shapefile;

			std::string temp_dir;

			bool overwrite_output;

		};

		extern Parameters* params;
	}
}
