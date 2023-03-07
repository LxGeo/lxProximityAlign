#include "parameters.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{
		Parameters::Parameters(int argc, char* argv[])
		{
			init();
			CLI::App app{ "lxProximityAlign" };
			app.add_option("--ishp", input_shapefile_to_align, "Polygons shapefile to align!")->required()->check(CLI::ExistingFile);
			app.add_option("--rshp", input_ref_shapefile, "Polygons shapefile used as reference!")->required()->check(CLI::ExistingFile);
			app.add_option("-o, --output", output_shapefile, "Output path of aligned polygons shapefile!");
			app.add_option("--r_imd", r_imd_path, "Metadata file respective rshp!")->required()->check(CLI::ExistingFile);
			app.add_option("--i_imd", i_imd_path, "Metadata file respective to ishp!")->required()->check(CLI::ExistingFile);
			app.add_flag("--keep_geometries", keep_geometries, "Flag to keep initial geometries!");
			app.add_option("--max_disparity", max_disparity, "Mximum search disparity in meteres!")->check(CLI::Range(1, 10000));

			try {
				\
					(app).parse((argc), (argv)); \
			}
			catch (const CLI::ParseError& e) {
				\
					(app).exit(e); \
			}
			
		}

		Parameters::~Parameters(){}

		bool Parameters::initialized()
		{
			return !input_shapefile_to_align.empty();
		}

		void Parameters::init()
		{
			printed_help = false;

			input_shapefile_to_align.clear();
			input_ref_shapefile.clear();
			r_imd_path.clear();
			i_imd_path.clear();

			output_basename = "result";
			output_shapefile = "result.shp";
			temp_dir = "/temp_dir/";
			overwrite_output = false;
			keep_geometries = false;
			max_disparity = 200;

		}

		Parameters* params = nullptr;
	}
}