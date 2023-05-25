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
			app.add_option("--couple_path", couple_path, "Json file containing epipolar creation parameters!")->check(CLI::ExistingFile);
			app.add_flag("--keep_geometries", keep_geometries, "Flag to keep initial geometries!");
			app.add_option("--max_disparity", max_disparity, "Mximum search disparity in meteres!")->check(CLI::Range(1, 10000));
			app.add_option("--ndbv", neighbour_distance_band_values, "A container of neighbour distance band values used to create components! \
				Ascending order distances, where the first value correspond to the last step distance (-1 value means disconnect all features)!");

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
			couple_path.clear();

			output_basename = "result";
			output_shapefile = "result.shp";
			temp_dir = "/temp_dir/";
			overwrite_output = false;
			keep_geometries = false;
			max_disparity = 150;
			neighbour_distance_band_values = { 0, 1, 10, 20, 40, 100 };

		}

		Parameters* params = nullptr;
	}
}