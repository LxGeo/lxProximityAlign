#include "optim_proximity_aligner.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		bool OptimProximityAligner::pre_check() {

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

			if (!params->couple_path.empty()) {
				std::ifstream f(params->couple_path);
				nlohmann::json data = nlohmann::json::parse(f);
				couple_rotation_angle = data["rotation_angle"];
				couple_v_displacement = data["v_disp"];
				std::vector<double> xy_values = data["xy_cst"]; 
				assert(xy_values.size() == 2 && "Epipolar couple file is corrupt!");
				xy_cst = std::make_pair(xy_values[0], xy_values[1]);
			}
			else if (!params->r_imd_path.empty()){
				IMetaData r_imd(params->r_imd_path);
				IMetaData i_imd(params->i_imd_path);
				xy_cst = compute_roof2roof_constants(RADS(i_imd.satAzimuth), RADS(i_imd.satElevation), RADS(r_imd.satAzimuth), RADS(r_imd.satElevation));
			}
			else {
				xy_cst = { 0,0 };
			}
			MAX_DISP = params->max_disparity;

			// Check neighbour_distance_band_values
			auto it = std::adjacent_find(params->neighbour_distance_band_values.begin(), params->neighbour_distance_band_values.end(), std::greater<double>());
			if (it != params->neighbour_distance_band_values.end())
			{
				std::cout << "Neighbour distance band values should be in acsending order!" << std::endl;
				return false;
			}

			//output dirs creation
			if (!endswith(params->output_shapefile, ".shp")) {
				std::cout << "Output shapefile should have .shp extension!" << std::endl;
				return false;
			}
			boost::filesystem::path output_path(params->output_shapefile);
			boost::filesystem::path output_parent_dirname = output_path.parent_path();
			boost::filesystem::path output_temp_path = output_parent_dirname / params->temp_dir;
			params->temp_dir = output_temp_path.string();

			boost::filesystem::create_directories(output_parent_dirname);
			boost::filesystem::create_directory(output_temp_path);

			if (boost::filesystem::exists(output_path) && !params->overwrite_output) {
				std::cout << "output shapefile already exists:!" << std::endl;
				std::cout << "Add --overwrite_output !" << std::endl;
				return false;
			}

			return true;

		}

		void OptimProximityAligner::run() {

			PolygonsShapfileIO target_shape, ref_shape;
			bool target_loaded = target_shape.load_shapefile(params->input_shapefile_to_align, true);
			bool ref_loaded = ref_shape.load_shapefile(params->input_ref_shapefile, true);
			// get common AOI extents
			OGREnvelope target_envelope, ref_envelope;
			target_shape.vector_layer->GetExtent(&target_envelope); ref_shape.vector_layer->GetExtent(&ref_envelope);
			OGREnvelope union_envelope(target_envelope); union_envelope.Merge(ref_envelope);
			union_envelope = box_buffer(union_envelope, std::abs(MAX_DISP));

			std::string out_proximity = (boost::filesystem::path(params->temp_dir) / "proximity.tif").string();
			polygons2proximity(ref_shape, out_proximity, &union_envelope, 0.5, 0.5, ProximityMapStrategy::contours);

			RasterIO ref_raster = RasterIO(out_proximity, GA_ReadOnly, false);
			GeoImage<cv::Mat> ref_gimg = GeoImage<cv::Mat>::from_file(out_proximity);
			float null_value = ref_gimg.no_data.value_or(FLT_MAX);
			RasterPixelsStitcher RPR(ref_gimg);

			std::unordered_map<std::string, matrix> matrices_map;
			matrices_map["proximity"] = ref_raster.raster_data;

			// Load shapefile
			Boost_Polygon_2 extents;
			//boost::geometry::read_wkt(
			//	"Polygon ((685997.13262533536180854 4630918.40684657730162144, 686254.43405542231630534 4630919.21343413181602955, 686252.8208803121233359 4630808.71093907952308655, 685999.149094223161228 4630811.13070174492895603, 685997.13262533536180854 4630918.40684657730162144))", extents);
			GeoVector<Boost_Polygon_2> in_gvector = GeoVector<Boost_Polygon_2>::from_file(params->input_shapefile_to_align, "", extents);
			OGRSpatialReference spatial_ref;
			VProfile vpr = VProfile::from_gdal_dataset(load_gdal_vector_dataset_shared_ptr(params->input_shapefile_to_align));
			spatial_ref.importFromWkt(vpr.s_crs_wkt.c_str());

			/*A function used in the optimizer to measure how well a geometry is aligned returning a value varying between 0-1 where 0 is not aligned and 1 is perfectly aligned
			*/
			
			auto fitness_from_stats_functor = [this](numcpp::DetailedStats<float>& stats)->float {
				if (stats.empty())
					return 1;

				auto sigmoid = [](const double& value, double offset = 0.0, double scale = 1.0, double exp_scale = 1.0) {return scale / (1.0 + std::exp(-exp_scale * value)) + offset; };

				const double TH = 3.0;
				auto countSubsequencesSum = [&TH](const std::vector<float>& values) {
					int count = 0;
					int subsequenceLength = 0;
					for (const auto& value : values) {
						if (value<TH) {
							++subsequenceLength;
						}
						else {
							count += (subsequenceLength > 1) ? subsequenceLength : 0;
							subsequenceLength = 0;
						}
					}
					count += (subsequenceLength > 1) ? subsequenceLength : 0;
					return count;
				};

				double alpha = 0.6, beta = 0.2, gamma = 0.2;
				auto& values_container = stats.values_vector();
				int under_th_count = countSubsequencesSum(values_container);
				double under_th_ratio = double(under_th_count) / values_container.size();
				double w_objective = alpha * sigmoid(stats.mean(), -1.0, 2.0, 0.1) + beta * (1 - under_th_ratio) + gamma * sigmoid(stats.stdev(), -1.0, 2.0, 0.1);
				return w_objective / (alpha+beta+gamma);
			};
			
			std::string DISP_COLUMN_NAME = "DISP";

			
			//nm_proximity_align_1d(matrices_map, ref_gimg, in_gvector, params->neighbour_distance_band_values, xy_cst, fitness_from_stats_functor, MAX_DISP, DISP_COLUMN_NAME);
			std::pair<std::string, std::string> disp_column_names = { "DISP_X", "DISP_Y" };
			pagmo_proximity_align_polygon(matrices_map, ref_gimg, in_gvector, params->neighbour_distance_band_values, fitness_from_stats_functor, MAX_DISP, disp_column_names);
			//nm_proximity_align(matrices_map, ref_gimg, in_gvector, params->neighbour_distance_band_values, fitness_from_stats_functor, MAX_DISP, disp_column_names);
			// Assign confidence and displacement			
			
			for (auto& c_gwa : in_gvector.geometries_container) {

				//double disp_value = c_gwa.get_double_attribute(DISP_COLUMN_NAME);
				//double ddx = disp_value * xy_cst.first;
				double ddx = c_gwa.get_double_attribute(disp_column_names.first); //disp_value * xy_cst.first;
				//double ddy = disp_value * xy_cst.second;
				double ddy = c_gwa.get_double_attribute(disp_column_names.second); //disp_value * xy_cst.second;
				bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(ddx, ddy);
				auto translated_geom = translate_geometry(c_gwa.get_definition(), trans_obj);

				if (!params->keep_geometries) {					
					c_gwa.set_definition(translated_geom);
				}

				auto  stitched_pixels = RPR.readPolygonPixels<float>(translated_geom, RasterPixelsStitcherStartegy::contours);

				std::vector<float> adj_diff; adj_diff.reserve(stitched_pixels.size());
				std::adjacent_difference(stitched_pixels.begin(), stitched_pixels.end(), std::back_inserter(adj_diff), [](float& a, float b) {return std::abs(a - b); });
				double volatility = (adj_diff.size()>1) ? std::accumulate(std::next(adj_diff.begin()), adj_diff.end(), 0.0) / (adj_diff.size()-1) : -1;
				auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
				c_gwa.set_double_attribute("confidence", fitness_from_stats_functor(stats));
				c_gwa.set_double_attribute("variance", stats.variance());
				c_gwa.set_double_attribute("coeff_var", stats.coeff_var());
				c_gwa.set_double_attribute("stdev", stats.stdev());
				c_gwa.set_double_attribute("mean", stats.mean());
				c_gwa.set_double_attribute("max", (!stats.empty()) ? stats.max(): DBL_MAX);
				c_gwa.set_double_attribute("min", (!stats.empty()) ? stats.min(): DBL_MAX);
				c_gwa.set_double_attribute("vola1", volatility);
			}

			std::cout << "Writing outfile!" << std::endl;
			in_gvector.to_file(params->output_shapefile, &spatial_ref);

		}

	}
}