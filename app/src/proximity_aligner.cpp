#include "proximity_aligner.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		bool ProximityAligner::pre_check() {

			std::cout << "Pre check input parameters" << std::endl;

			PolygonsShapfileIO target_shape, ref_shape;
			bool target_loaded = target_shape.load_shapefile(params->input_shapefile_to_align, true);
			if (!target_loaded) {
				std::cout << "Error loading shapefile to align at: " << params->input_shapefile_to_align << std::endl;
				return false;
			}
			bool ref_loaded = ref_shape.load_shapefile(params->input_ref_shapefile, true);
			if (!ref_loaded) {
				std::cout << "Error loading reference shapefile at: " << params->input_ref_shapefile << std::endl;
				return false;
			}

			if (!target_shape.spatial_refrence->IsSame(ref_shape.spatial_refrence)) {
				std::cout << "Input shapefiles have different spatial reference system!" << std::endl;
				return false;
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
				boost::filesystem::create_directory(output_temp_path);
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

			PolygonsShapfileIO target_shape, ref_shape;
			bool target_loaded = target_shape.load_shapefile(params->input_shapefile_to_align, true);
			bool ref_loaded = ref_shape.load_shapefile(params->input_ref_shapefile, true);
			// get common AOI extents
			OGREnvelope target_envelope, ref_envelope;
			target_shape.vector_layer->GetExtent(&target_envelope); ref_shape.vector_layer->GetExtent(&ref_envelope);
			OGREnvelope union_envelope(target_envelope); union_envelope.Merge(ref_envelope);
						
			std::string out_proximity = (boost::filesystem::path(params->temp_dir) / "proximity.tif").string();
			polygons2proximity(ref_shape, out_proximity, &union_envelope, 0.5,0.5,ProximityMapStrategy::contours);

			RasterIO ref_raster = RasterIO(out_proximity, GA_ReadOnly, false);

			std::map<std::string, matrix> matrices_map;
			matrices_map["proximity"] = ref_raster.raster_data;
			matrix grad_x; cv::Sobel(ref_raster.raster_data, grad_x, CV_32FC1, 0, 1);
			matrix grad_y; cv::Sobel(ref_raster.raster_data, grad_y, CV_32FC1, 1, 0);
			if (true) { //save gradients
				RasterIO grad_x_raster = RasterIO(ref_raster, grad_x);
				RasterIO grad_y_raster = RasterIO(ref_raster, grad_y);
				std::string grad_x_out_path = (boost::filesystem::path(params->temp_dir) / "grad_x.tif").string();
				std::string grad_y_out_path = (boost::filesystem::path(params->temp_dir) / "grad_y.tif").string();
				grad_x_raster.write_raster(grad_x_out_path, false);
				grad_y_raster.write_raster(grad_y_out_path, false);

			}
			matrices_map["grad_x"] = grad_x; //RasterIO(params->x_g_raster_path, GA_ReadOnly, false);
			matrices_map["grad_y"] = grad_y; //RasterIO(params->y_g_raster_path, GA_ReadOnly, false);

			
			// Load shapefile
			GeoVector<Boost_Polygon_2> in_gvector = GeoVector<Boost_Polygon_2>::from_file(params->input_shapefile_to_align);
			OGRSpatialReference spatial_ref;
			VProfile vpr = VProfile::from_gdal_dataset(load_gdal_vector_dataset_shared_ptr(params->input_shapefile_to_align));
			spatial_ref.importFromWkt(vpr.s_crs_wkt.c_str());
			
			std::vector<Boost_Polygon_2> aligned_polygon= alignmentNeighbour(matrices_map, ref_raster, in_gvector.geometries_container);


			PolygonsShapfileIO aligned_out_shapefile = PolygonsShapfileIO(params->output_shapefile, &spatial_ref);
			auto polygons_with_attrs = transform_to_geom_with_attr<Boost_Polygon_2>(aligned_polygon);
			std::cout << "Writing outfile!" << std::endl;
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