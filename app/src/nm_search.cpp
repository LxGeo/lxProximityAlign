#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "nm_search.h"
#include "optim.hpp"
#include "numcpp/line_space.h"

#include "graph_weights/spatial_weights.h"
#include "affine_geometry/affine_transformer.h"
//#include "raster_stitch.h"
#include "stitching/vector_on_raster_stitcher.h"
#include "io_shapefile.h"
#include "parameters.h"
#include "numcpp/stats.h"


namespace LxGeo
{
	namespace lxProximityAlign
	{

		void nm_proximity_align(
			std::unordered_map<std::string, matrix>& matrices_map, GeoImage<cv::Mat>& ref_gimg,
			GeoVector<Boost_Polygon_2>& input_geovector,
			std::pair<std::string, std::string> OBJECTIVE_FIELD_NAME_PAIR) {

			std::vector<double> neighbour_distance_band_values = {0.01, 1, 5, 10, 20};//numcpp::linspace(0.0, 10, 9);
			SpatialWeights<Boost_Polygon_2> PSW = SpatialWeights<Boost_Polygon_2>::from_geovector(input_geovector);
			WeightsDistanceBandParams wdbp = { 20, false, -1, [](double x)->double { return x; } };
			PSW.fill_distance_band_graph(wdbp);

			RasterPixelsStitcher RPR(ref_gimg);
			float null_value = ref_gimg.no_data.value_or(FLT_MAX);

			auto objective_fn = [&RPR, &null_value, &OBJECTIVE_FIELD_NAME_PAIR](const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)->double {
				
				double obj_val = 0;
				double OUT_OF_BOUNDS_PENALTY = 100.0;

				auto resp_gwas = static_cast<std::list<Geometries_with_attributes<Boost_Polygon_2>*>*>(opt_data);

				for (auto& c_gwa : *resp_gwas) {

					double c_geom_disp_x = vals_inp(0) + c_gwa->get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first);
					double c_geom_disp_y = vals_inp(1) + c_gwa->get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second);

					bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(c_geom_disp_x, c_geom_disp_y);
					auto translated_geometry = translate_geometry(c_gwa->get_definition(), trans_obj);
					auto  stitched_pixels = RPR.readPolygonPixels<float>(translated_geometry, RasterPixelsStitcherStartegy::contours);
					auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
					double norm_sum;
					if (!stats.empty())
						norm_sum = stats.sum() / stats.count();
					else
						norm_sum = OUT_OF_BOUNDS_PENALTY;
					obj_val += norm_sum;
				}
				return obj_val;
				
			};

			
			double MAX_DISP=params->max_disparity;
			for (auto distance_val_iter = neighbour_distance_band_values.rbegin(); distance_val_iter != neighbour_distance_band_values.rend(); ++distance_val_iter) {

				MAX_DISP /= 2;

				auto disconnection_lambda = [&distance_val_iter](double x)->bool {return x>*distance_val_iter;};
				PSW.disconnect_edges(disconnection_lambda);
				PSW.run_labeling();
				std::cout << "Aligning N components: " << PSW.n_components << std::endl;

				std::map<size_t, std::list<size_t>> components_polygons_map;
				std::map<size_t, Eigen::VectorXd> components_transforms_map;
				for (size_t _idx = 0; _idx < PSW.component_labels.size(); ++_idx) {
					size_t comp_id = PSW.component_labels[_idx];
					if (components_polygons_map.find(comp_id) != components_polygons_map.end())
						components_polygons_map[comp_id].push_back(_idx);
					else
					{
						components_polygons_map[comp_id] = std::list<size_t>();
						components_polygons_map[comp_id].push_back(_idx);
						components_transforms_map[comp_id] = Eigen::VectorXd(2); 
						components_transforms_map[comp_id] << MAX_DISP * (double)rand() / RAND_MAX, MAX_DISP* (double)rand() / RAND_MAX;// 0.5,0.4;
					}
				}

				tqdm bar;
				for (size_t comp_idx = 0; comp_idx < PSW.n_components; ++comp_idx) {
					bar.progress(comp_idx, PSW.n_components);
					auto& init_vals = components_transforms_map[comp_idx];
					std::list<Geometries_with_attributes<Boost_Polygon_2>*> respective_polygons;
					for (auto poly_idx : components_polygons_map[comp_idx]) respective_polygons.push_back(&input_geovector.geometries_container[poly_idx]);
					
					optim::algo_settings_t c_settings;
					c_settings.iter_max = 50;
					c_settings.rel_objfn_change_tol = 5;
					c_settings.vals_bound = true; c_settings.lower_bounds = -MAX_DISP * Eigen::VectorXd::Ones(2); c_settings.upper_bounds = MAX_DISP * Eigen::VectorXd::Ones(2);
					//c_settings.print_level = 3;
					optim::Mat_t simplex_points(3, 2);
					simplex_points.row(0) << MAX_DISP * ((double)rand() / RAND_MAX - 0.5), MAX_DISP* ((double)rand() / RAND_MAX - 0.5); // 1.5, 1.5;
					simplex_points.row(1) << MAX_DISP * ((double)rand() / RAND_MAX - 0.5), MAX_DISP* ((double)rand() / RAND_MAX - 0.5); //0.1, 0.15;
					simplex_points.row(2) << MAX_DISP * ((double)rand() / RAND_MAX - 0.5), MAX_DISP* ((double)rand() / RAND_MAX - 0.5); //-1.0, 1.0;
					c_settings.nm_settings.custom_initial_simplex = true; c_settings.nm_settings.initial_simplex_points = simplex_points;

					bool success = optim::nm(init_vals, objective_fn, &respective_polygons, c_settings);
					if (!success) {
						std::cout << "optimization failed!" << std::endl; continue;
					}
					else {
						for (auto& c_gwa : respective_polygons) {
							c_gwa->set_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first, init_vals(0) + c_gwa->get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first));
							c_gwa->set_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second, init_vals(1) + c_gwa->get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second));
						}
					}
				}
				bar.finish();
			}

		}

		/*
		Given a geovector, a proximity map and epipolar direction, assign disparity value to each polygon corresponding to the best alignment.
		Geometries are intact.
		*/

		void nm_proximity_align_1d(
			std::unordered_map<std::string, matrix>& matrices_map, GeoImage<cv::Mat>& ref_gimg,
			GeoVector<Boost_Polygon_2>& input_geovector,
			std::pair<double, double>& r2r_constants,
			std::string OBJECTIVE_FIELD_NAME
		) {


			double MAX_DISP = params->max_disparity;
			
			std::vector<double> neighbour_distance_band_values = { 1, 5, 10, 20 };
			SpatialWeights<Boost_Polygon_2> PSW = SpatialWeights<Boost_Polygon_2>::from_geovector(input_geovector);
			WeightsDistanceBandParams wdbp = { neighbour_distance_band_values[neighbour_distance_band_values.size()-1], false, -1, [](double x)->double { return x; } };
			PSW.fill_distance_band_graph(wdbp);

			RasterPixelsStitcher RPR(ref_gimg);
			float null_value = ref_gimg.no_data.value_or(FLT_MAX);

			auto geometry_alignment_evaluator = [&RPR, &null_value](const Boost_Polygon_2& ref_geometry)->double {
				auto stitched_pixels = RPR.readPolygonPixels<float>(ref_geometry, RasterPixelsStitcherStartegy::contours);

				// Measure alignment continuity where this value represents the percentage of contiguous pixels that are under a defined threshhold TH (Higher -> well aligned)
				// Measure alignment completness where this value correspond to the percentage of values that are under a defined threshhold TH (Higher -> well aligned)
				int consecutive_touched_sum = 0, consecutive_touched_max = 0;
				int c_consecutive = 0;
				double TH = 5.0;
				for (const auto& c_value : stitched_pixels) {
					if (c_value < TH) c_consecutive += 1;
					else {
						consecutive_touched_sum += c_consecutive;
						consecutive_touched_max = (c_consecutive > consecutive_touched_max) ? c_consecutive : consecutive_touched_max;
						c_consecutive = 0;
					}
				}
				consecutive_touched_sum += c_consecutive;
				consecutive_touched_max = (c_consecutive > consecutive_touched_max) ? c_consecutive : consecutive_touched_max;
				double alignment_continuity = (stitched_pixels.empty()) ? -1.0 : double(consecutive_touched_max) / stitched_pixels.size();
				double alignment_completnes = (stitched_pixels.empty()) ? -1.0 : double(consecutive_touched_sum) / stitched_pixels.size();

				// Measure detailed stats
				auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);

				// The combined objective value
				double objective = (stats.empty()) ? DBL_MAX : stats.mean(); //* stats.variance() / alignment_continuity / alignment_completnes;
				return objective;
			};

			// init disparity
			for (auto& c_gwa : input_geovector.geometries_container) {
				c_gwa.set_double_attribute(OBJECTIVE_FIELD_NAME, 0.0);
				c_gwa.set_double_attribute("OBJ", geometry_alignment_evaluator(c_gwa.get_definition()));
			}

			auto objective_fn = [&RPR, &r2r_constants, &null_value, &OBJECTIVE_FIELD_NAME, &geometry_alignment_evaluator](const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)->double {
				double obj_val = 0;
				double OUT_OF_BOUNDS_PENALTY = 1e5;				

				auto resp_gwas = static_cast<std::list<Geometries_with_attributes<Boost_Polygon_2>*>*>(opt_data);

				auto c_component_height = vals_inp(0);				
					
				for (auto& c_gwa : *resp_gwas) {

					double c_geom_height = c_gwa->get_double_attribute(OBJECTIVE_FIELD_NAME) + c_component_height;
					bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(r2r_constants.first * c_geom_height, r2r_constants.second * c_geom_height);
					auto translated_geometry = translate_geometry(c_gwa->get_definition(), trans_obj);
					double c_geometry_obj_val = geometry_alignment_evaluator(translated_geometry);
					
					if (c_gwa->get_double_attribute("OBJ") > c_geometry_obj_val) {
						c_gwa->set_double_attribute("OBJ", c_geometry_obj_val);
						c_gwa->set_double_attribute(OBJECTIVE_FIELD_NAME, c_geom_height);
					}					
					obj_val += c_geometry_obj_val;
				}
				return obj_val;
			};
			
			for (auto distance_val_iter = neighbour_distance_band_values.rbegin(); distance_val_iter != neighbour_distance_band_values.rend(); ++distance_val_iter) {

				MAX_DISP /= 2;

				auto disconnection_lambda = [&distance_val_iter](double x)->bool {return x > *distance_val_iter; };
				PSW.disconnect_edges(disconnection_lambda);
				PSW.run_labeling();
				std::cout << "Aligning N components: " << PSW.n_components << std::endl;

				std::map<size_t, std::list<size_t>> components_polygons_map;
				std::map<size_t, Eigen::VectorXd> components_transforms_map;
				for (size_t _idx = 0; _idx < PSW.component_labels.size(); ++_idx) {
					size_t comp_id = PSW.component_labels[_idx];
					if (components_polygons_map.find(comp_id) != components_polygons_map.end())
						components_polygons_map[comp_id].push_back(_idx);
					else
					{
						components_polygons_map[comp_id] = std::list<size_t>();
						components_polygons_map[comp_id].push_back(_idx);
						components_transforms_map[comp_id] = Eigen::VectorXd(1);
						components_transforms_map[comp_id] << MAX_DISP * (double)rand() / RAND_MAX;
					}
				}

				tqdm bar;
				GeoVector<Boost_Polygon_2> component_gvector;
				for (size_t comp_idx = 0; comp_idx < PSW.n_components; ++comp_idx) {
					bar.progress(comp_idx, PSW.n_components);
					auto& init_vals = components_transforms_map[comp_idx];
					std::list<Geometries_with_attributes<Boost_Polygon_2>*> respective_polygons;
					for (auto poly_idx : components_polygons_map[comp_idx]) respective_polygons.push_back(&input_geovector.geometries_container[poly_idx]);

					// Init optimizer parameters
					optim::algo_settings_t c_settings;
					c_settings.iter_max = 50;
					c_settings.rel_objfn_change_tol = 5;
					c_settings.vals_bound = true; c_settings.lower_bounds = -3 * Eigen::VectorXd::Ones(1); c_settings.upper_bounds = MAX_DISP * Eigen::VectorXd::Ones(1);
					optim::Mat_t simplex_points(2, 1);
					simplex_points.row(0) << MAX_DISP * ((double)rand() / RAND_MAX ); 
					simplex_points.row(1) << MAX_DISP * ((double)rand() / RAND_MAX ); 
					c_settings.nm_settings.custom_initial_simplex = true; c_settings.nm_settings.initial_simplex_points = simplex_points;

					bool success = optim::nm(init_vals, objective_fn, &respective_polygons, c_settings);
					if (!success) { 
						std::cout << "optimization failed!" << std::endl; continue;
					}
					else {
						for (auto& c_gwa : respective_polygons) {
							//c_gwa->set_double_attribute(OBJECTIVE_FIELD_NAME, init_vals(0) + c_gwa->get_double_attribute(OBJECTIVE_FIELD_NAME));
							component_gvector.add_geometry(*c_gwa);
							auto saved_gwa = component_gvector.geometries_container.rbegin();
							saved_gwa->set_int_attribute("cmp_id", comp_idx);
							double disp_value = saved_gwa->get_double_attribute(OBJECTIVE_FIELD_NAME);
							double ddx = disp_value * r2r_constants.first;
							double ddy = disp_value * r2r_constants.second;
							bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(ddx, ddy);
							//saved_gwa->set_definition(translate_geometry(saved_gwa->get_definition(), trans_obj));
						}
					}
				}
				bar.finish();
				std::string comp_file = params->temp_dir + "//comp_" + std::to_string(*distance_val_iter) + ".shp";
				auto geometry_shifter_functor = [&OBJECTIVE_FIELD_NAME, &r2r_constants](const Geometries_with_attributes<Boost_Polygon_2>& c_gwa) {
					double disp_value = c_gwa.get_double_attribute(OBJECTIVE_FIELD_NAME);
					double ddx = disp_value * r2r_constants.first;
					double ddy = disp_value * r2r_constants.second;
					bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(ddx, ddy);
					auto translated_geometry = translate_geometry(c_gwa.get_definition(), trans_obj);
					auto new_gwa = Geometries_with_attributes<Boost_Polygon_2>(c_gwa);
					new_gwa.set_definition(translated_geometry);
					return new_gwa;
				};
				
				GeoVector<Boost_Polygon_2> to_save_components; transform_geovector<Boost_Polygon_2, Boost_Polygon_2>(component_gvector, to_save_components, geometry_shifter_functor);
				to_save_components.to_file(comp_file);
			}
		}

	}
}