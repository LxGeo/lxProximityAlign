#include "nm_search.h"

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
			std::vector<double>& neighbour_distance_band_values,
			std::pair<double, double>& r2r_constants,
			std::function<float(numcpp::DetailedStats<float>&)> fitness_from_stats_functor,
			double MAX_DISP,
			std::string DISP_COLUMN_NAME
		) {
			std::string FITT_COLUMN_NAME = "FITT";
			double FITNESS_THRESH = 0.0;
			// Init DISP field wiith zero values
			for (auto& c_gwa : input_geovector.geometries_container) {
				c_gwa.set_double_attribute(DISP_COLUMN_NAME, 0);
				c_gwa.set_double_attribute(FITT_COLUMN_NAME, FLT_MAX);
			}
			
			// Creating spatial weights from input geometries
			SpatialWeights<Boost_Polygon_2> PSW = SpatialWeights<Boost_Polygon_2>::from_geovector(input_geovector);
			WeightsDistanceBandParams wdbp = { neighbour_distance_band_values[neighbour_distance_band_values.size()-1], false, -1, [](double x)->double { return x; } };
			PSW.fill_distance_band_graph(wdbp);

			// Loading proximity map stitcher
			RasterPixelsStitcher RPR(ref_gimg);
			float null_value = ref_gimg.no_data.value_or(FLT_MAX);

			/*****
			*Lambda function used to return the fitness of a geometry within the proximity map 
			*****/
			auto geometry_fitness_evaluator = [&RPR, &null_value, &fitness_from_stats_functor](const Boost_LineString_2& ref_geometry)->double {
				auto stitched_pixels = RPR.readLineStringPixels<float>(ref_geometry, RasterPixelsStitcherStartegy::contours);
				// Measure detailed stats
				auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
				return fitness_from_stats_functor(stats);
			};
			
			int pyramid_level = neighbour_distance_band_values.size();
			for (auto distance_val_iter = neighbour_distance_band_values.rbegin(); distance_val_iter != neighbour_distance_band_values.rend(); ++distance_val_iter, pyramid_level--) {

				auto disconnection_lambda = [&distance_val_iter](double x)->bool {return x - 1e-3 > *distance_val_iter; };
				PSW.disconnect_edges(disconnection_lambda);
				PSW.run_labeling();
				std::cout << "Aligning " << PSW.n_components << " components! " << std::endl;

				std::map<size_t, std::list<size_t>> components_polygons_map;
				for (size_t _idx = 0; _idx < PSW.component_labels.size(); ++_idx) {
					double c_fitness_value = input_geovector.geometries_container[_idx].get_double_attribute(FITT_COLUMN_NAME);
					if (c_fitness_value < FITNESS_THRESH)
						continue;

					size_t comp_id = PSW.component_labels[_idx];
					if (components_polygons_map.find(comp_id) != components_polygons_map.end())
						components_polygons_map[comp_id].push_back(_idx);
					else
					{
						components_polygons_map[comp_id] = std::list<size_t>();
						components_polygons_map[comp_id].push_back(_idx);
					}
				}

				tqdm bar;
				GeoVector<Boost_Polygon_2> component_gvector;
				// Aligning components one by one
				for (size_t comp_idx = 0; comp_idx < PSW.n_components; ++comp_idx) {
					if (components_polygons_map[comp_idx].empty())
						continue;
					// Get previous disparity value assigned to the broader component
					double previous_step_disp_value = input_geovector.geometries_container[*components_polygons_map[comp_idx].begin()].get_double_attribute(DISP_COLUMN_NAME);
					bar.progress(comp_idx, PSW.n_components);
					std::list<Geometries_with_attributes<Boost_Polygon_2>*> respective_polygons;
					for (auto poly_idx : components_polygons_map[comp_idx]) {
						respective_polygons.push_back(&input_geovector.geometries_container[poly_idx]);
					}
					Arrangement respective_polygons_arrangment = ArrangmentFromPolygons(respective_polygons);
					
					// Extract component exterior edges using arrangment
					std::list<Boost_LineString_2> component_exterior_edges;
					
					/*for (auto eit = respective_polygons_arrangment.edges_begin(); eit != respective_polygons_arrangment.edges_end(); ++eit) {
						auto& src = eit->curve().source();
						auto& tar = eit->curve().target();
						auto c_segment = Boost_LineString_2({ { CGAL::to_double(src.x()), CGAL::to_double(src.y()) }, { CGAL::to_double(tar.x()), CGAL::to_double(tar.y()) } });
						component_exterior_edges.push_back(c_segment);
					};*/
					
					for (auto eit = respective_polygons_arrangment.halfedges_begin(); eit != respective_polygons_arrangment.halfedges_end(); ++eit) {
						// ignore if edge is shared
						if (!eit->twin()->face()->is_unbounded())
							continue;
						auto& src = eit->curve().source();
						auto& tar = eit->curve().target();
						component_exterior_edges.push_back(
							Boost_LineString_2({ { CGAL::to_double(src.x()), CGAL::to_double(src.y()) }, { CGAL::to_double(tar.x()), CGAL::to_double(tar.y()) } })
						);
					};
					

					auto bound_objective = [&RPR, &r2r_constants, &null_value, &geometry_fitness_evaluator, &component_exterior_edges](double c_component_height)->double {
						double total_fitness = 0;
						double total_length = 0;
						for (auto& c_exterior_edge : component_exterior_edges) {
							double c_edge_length = bg::length(c_exterior_edge);
							bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(r2r_constants.first * c_component_height, r2r_constants.second * c_component_height);
							auto translated_geometry = translate_geometry(c_exterior_edge, trans_obj);
							double c_geometry_obj_val = geometry_fitness_evaluator(translated_geometry);
							total_fitness += (c_geometry_obj_val* c_edge_length);
							total_length += c_edge_length;
						}
						double mean_fitness = total_fitness / total_length;						
						return mean_fitness;
					};;
					auto cached_objective = FunctionCache<double, double>(bound_objective);

					double min_search_bound = - MAX_DISP/ pyramid_level /2 + previous_step_disp_value;
					double max_search_bound =  MAX_DISP + previous_step_disp_value;
					double optimal_disparity = custom_splitting_search(cached_objective, min_search_bound, max_search_bound, 0.1, 500, (max_search_bound- min_search_bound)/10);
					//double optimal_disparity = particleSwarmOptimization(cached_objective, min_search_bound, max_search_bound, 10,50);
					/*double optimal_disparity = simulated_annealing_bounded(
						cached_objective, min_search_bound, max_search_bound,
						(max_search_bound + min_search_bound) / 2,
						1.0, 0.95, 100);*/
					double optimal_fitness = cached_objective(optimal_disparity);
					double previous_optimal_fitness = cached_objective(0.0);

					// ignore single geometry alignment if optimal_fitness is over a threshhold
					double single_geom_fitt_threshhold = 0.1;
					if (components_polygons_map[comp_idx].size() == 1 && optimal_fitness > single_geom_fitt_threshhold) {
						optimal_disparity = previous_step_disp_value;
						optimal_fitness = previous_optimal_fitness;
					}

					// Ignore updating if fitness is high
					if (optimal_fitness > previous_optimal_fitness || optimal_disparity<0) {
						//std::cout << "optimization failed!" << std::endl; continue;
						optimal_disparity = previous_step_disp_value;
						optimal_fitness = previous_optimal_fitness;
					}

					
					for (auto c_poly_ptr : respective_polygons) {
						
						double ddx = optimal_disparity * r2r_constants.first;
						double ddy = optimal_disparity * r2r_constants.second;
						bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(ddx, ddy);
						auto translated_geometry = translate_geometry(c_poly_ptr->get_definition(), trans_obj);
						auto stitched_pixels = RPR.readPolygonPixels<float>(translated_geometry, RasterPixelsStitcherStartegy::contours);
						auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
						double c_polygon_fitness = fitness_from_stats_functor(stats);
						double c_polygon_previous_fitness = c_poly_ptr->get_double_attribute(FITT_COLUMN_NAME);
						double c_polygon_previous_disparity = c_poly_ptr->get_double_attribute(DISP_COLUMN_NAME);

						bool alignment_success = c_polygon_fitness < c_polygon_previous_fitness;

						// add to component temporary shp
						if (!alignment_success) {
							bg::strategy::transform::translate_transformer<double, 2, 2> prev_trans_obj(
								c_polygon_previous_disparity * r2r_constants.first,
								c_polygon_previous_disparity* r2r_constants.second
							);
							component_gvector.add_geometry(translate_geometry(c_poly_ptr->get_definition(), prev_trans_obj));
						}
						else
							component_gvector.add_geometry( translated_geometry);

						auto saved_gwa = component_gvector.geometries_container.rbegin();
						saved_gwa->set_double_attribute("obj", optimal_fitness);
						saved_gwa->set_int_attribute("cmp_id", comp_idx);
						saved_gwa->set_double_attribute("disp", (alignment_success) ? optimal_disparity: c_polygon_previous_disparity);					
						saved_gwa->set_double_attribute("ind_obj", (alignment_success) ? c_polygon_fitness : c_polygon_previous_fitness);

						// updating fields
						c_poly_ptr->set_double_attribute(DISP_COLUMN_NAME, (alignment_success) ? optimal_disparity : c_polygon_previous_disparity);
						c_poly_ptr->set_double_attribute(FITT_COLUMN_NAME, (alignment_success) ? c_polygon_fitness : c_polygon_previous_fitness);
					}
					
				}
				bar.finish();
				std::string comp_file = params->temp_dir + "//comp_" + std::to_string(*distance_val_iter) + ".shp";				
				component_gvector.to_file(comp_file);
			}
		}

	}
}