#include "io_raster.h"
#include "io_shapefile.h"
#include "defs.h"
#include "parameters.h"
#include "proximity_aligner.h"
#include "raster_stitch.h"
#include "graph_weights/polygons_spatial_weights.h"
#include "graph_weights/spatial_weights.h"
#include "geometries_with_attributes/linestring_with_attributes.h"
#include "alignment_basic_optim.h"
#include "alignment_iterative_support_optim.h"
#include "alignment_neighbour_optim.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		bool ProximityAligner::pre_check() {

			std::cout << "Pre check input parameters" << std::endl;

			RasterIO ref_raster = RasterIO();
			if (!ref_raster.load_raster(params->input_proximity_raster_path.c_str())) {
				std::cout << "Cannot load proximity raster at: " << params->input_proximity_raster_path.c_str() << std::endl;
				return false;
			}
			
			std::list<std::string> rasters_list = { params->input_proximity_raster_path, params->x_g_raster_path, params->y_g_raster_path };
			for (auto& c_raster_path : rasters_list){
				RasterIO temp_raster = RasterIO();
				// correct raster path check
				if (!temp_raster.load_raster(c_raster_path.c_str())) {
					std::cout << "Cannot load raster at: " << c_raster_path.c_str() << std::endl;
					return false;
				}
				if (!temp_raster.compare(ref_raster)) {
					std::cout << "Raster at: " << c_raster_path.c_str() << " doesn't match proximity raster!" << std::endl;
					return false;
				}


			}

			//output dirs creation
			boost::filesystem::path output_path(params->output_basename);
			boost::filesystem::path output_parent_dirname = output_path.parent_path();
			boost::filesystem::path output_temp_path = output_parent_dirname / params->temp_dir;
			params->temp_dir = output_temp_path.string();
			if (!boost::filesystem::exists(output_parent_dirname))
			{
				boost::filesystem::create_directory(output_parent_dirname);
				std::cout << fmt::format("Directory Created: {}", output_parent_dirname.string()) << std::endl;
			}

			if (!boost::filesystem::exists(output_temp_path))
			{
				std::cout << fmt::format("Directory Created: {}", output_temp_path.string()) << std::endl;
			}

			if (boost::filesystem::exists(output_path) && !params->overwrite_output) {
				std::cout << fmt::format("output shapefile already exists: {}!", output_path.string()) << std::endl;
				std::cout << fmt::format("Add --overwrite_output !", output_path.string()) << std::endl;
				return false;
			}

			return true;

		}

		void ProximityAligner::run() {

			std::map<std::string, RasterIO> rasters_map;
			rasters_map["proximity"] = RasterIO(params->input_proximity_raster_path, GA_ReadOnly, false);
			rasters_map["grad_x"] = RasterIO(params->x_g_raster_path, GA_ReadOnly, false);
			rasters_map["grad_y"] = RasterIO(params->y_g_raster_path, GA_ReadOnly, false);
						
			RasterIO& ref_raster = rasters_map["proximity"];
			
			// Load shapefile
			PolygonsShapfileIO sample_shape;
			bool loaded = sample_shape.load_shapefile(params->input_shapefile_to_align, false);
			
			std::vector<Boost_Polygon_2> aligned_polygon= alignmentNeighbour(rasters_map, sample_shape.geometries_container);


			PolygonsShapfileIO aligned_out_shapefile = PolygonsShapfileIO(params->output_shapefile, sample_shape.spatial_refrence);
			auto polygons_with_attrs = transform_to_geom_with_attr<Boost_Polygon_2>(aligned_polygon);
			aligned_out_shapefile.write_shapefile(polygons_with_attrs);

			/*PolygonSpatialWeights PSW = PolygonSpatialWeights(sample_shape.geometries_container);
			WeightsDistanceBandParams wdbp = { 1, false, -1, [](double x)->double { return 1.0 / (1.0 + x); } };
			PSW.fill_distance_band_graph(wdbp);
			PSW.run_labeling();
			std::vector<LineString_with_attributes> edges_line_strings = PSW.export_edge_graph_as_LSwithAttr();
			LineStringShapfileIO edges_out_shapefile = LineStringShapfileIO(params->output_shapefile, sample_shape.spatial_refrence);
			edges_out_shapefile.write_linestring_shapefile(edges_line_strings);
			bool a = true;*/
		}

	}
}