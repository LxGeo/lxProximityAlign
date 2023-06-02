#include "nm_search.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{

		void nm_proximity_align(
			std::unordered_map<std::string, matrix>& matrices_map, GeoImage<cv::Mat>& ref_gimg,
			GeoVector<Boost_Polygon_2>& input_geovector,
			std::vector<double>& neighbour_distance_band_values,
			std::function<float(numcpp::DetailedStats<float>&)> fitness_from_stats_functor,
			double MAX_DISP,
			std::pair<std::string, std::string> OBJECTIVE_FIELD_NAME_PAIR) {

			std::string FITT_COLUMN_NAME = "FITT";

			for (auto& c_gwa : input_geovector.geometries_container) {
				c_gwa.set_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first, 0);
				c_gwa.set_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second, 0);
				c_gwa.set_double_attribute(FITT_COLUMN_NAME, 1);
			}
			SpatialWeights<Boost_Polygon_2> PSW = SpatialWeights<Boost_Polygon_2>::from_geovector(input_geovector);
			WeightsDistanceBandParams wdbp = { neighbour_distance_band_values[neighbour_distance_band_values.size() - 1], false, -1, [](double x)->double { return x; } };
			PSW.fill_distance_band_graph(wdbp);

			RasterPixelsStitcher RPR(ref_gimg);
			float null_value = ref_gimg.no_data.value_or(FLT_MAX);

			/*****
			*Lambda function used to return the fitness of a geometry within the proximity map
			*****/
			auto linestring_fitness_evaluator = [&RPR, &null_value, &fitness_from_stats_functor](const Boost_LineString_2& ref_geometry)->double {
				auto stitched_pixels = RPR.readLineStringPixels<float>(ref_geometry, RasterPixelsStitcherStartegy::contours);
				// Measure detailed stats
				auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
				return fitness_from_stats_functor(stats);
			};

			ND::PT major_vector = { 0.0,0.0 };
			
			for (auto distance_val_iter = neighbour_distance_band_values.rbegin(); distance_val_iter != neighbour_distance_band_values.rend(); ++distance_val_iter) {
				GeoVector<Boost_Polygon_2> component_gvector;
				//MAX_DISP /= 2;

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
					}
				}

				tqdm bar;
				for (size_t comp_idx = 0; comp_idx < PSW.n_components; ++comp_idx) {
					bar.progress(comp_idx, PSW.n_components);
					Eigen::VectorXd optimal_displacement(2); optimal_displacement << 0.0, 0.0;
					
					double previous_step_disp_x_value = input_geovector.geometries_container[*components_polygons_map[comp_idx].begin()].get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first);
					double previous_step_disp_y_value = input_geovector.geometries_container[*components_polygons_map[comp_idx].begin()].get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second);
					Eigen::VectorXd previous_step_disp_xy(2); previous_step_disp_xy << previous_step_disp_x_value, previous_step_disp_y_value;

					std::list<Geometries_with_attributes<Boost_Polygon_2>*> respective_polygons;
					for (auto poly_idx : components_polygons_map[comp_idx]) {
						respective_polygons.push_back(&input_geovector.geometries_container[poly_idx]);
					}
					Arrangement respective_polygons_arrangment = ArrangmentFromPolygons(respective_polygons);
					// Extract component exterior edges using arrangment
					std::list<Boost_LineString_2> component_exterior_edges;
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

					std::function<double(const std::pair<double,double>)> bound_objective = [&RPR, &null_value, &linestring_fitness_evaluator, &component_exterior_edges](const std::pair<double, double>& disp_xy)->double {
						double total_fitness = 0;
						double total_length = 0;
						for (auto& c_exterior_edge : component_exterior_edges) {
							double c_edge_length = bg::length(c_exterior_edge);
							bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(disp_xy.first, disp_xy.second);
							auto translated_geometry = translate_geometry(c_exterior_edge, trans_obj);
							double c_geometry_obj_val = linestring_fitness_evaluator(translated_geometry);
							total_fitness += (c_geometry_obj_val * c_edge_length);
							total_length += c_edge_length;
						}
						double mean_fitness = total_fitness / total_length;
						return mean_fitness;
					};;
					struct PairHash {						
						std::size_t operator () (const std::pair<double, double>& pair) const {
							std::size_t h1 = std::hash<double>()(pair.first);
							std::size_t h2 = std::hash<double>()(pair.second);
							return h1 ^ h2;
						}
					};
					auto cached_objective = FunctionCache<std::pair<double, double>, double, PairHash>(bound_objective);
					auto optim_adapted_objective_fn = [&cached_objective](const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)->double {
						return cached_objective({vals_inp(0), vals_inp(1)});
					};
					

					//bool success = optim::nm(optimal_displacement, optim_adapted_objective_fn, &respective_polygons, c_settings);
					std::function<double(ND::PT)> ND_adaptetd_f = [&cached_objective, &major_vector, &MAX_DISP](ND::PT pt) {

						auto cosine_similarity = [](const ND::PT& v1, const ND::PT& v2) {
							double dot = 0.0, denom_a = 0.0, denom_b = 0.0;

							for (unsigned int i = 0; i < v1.size(); ++i) {
								dot += v1[i] * v2[i];
								denom_a += v1[i] * v1[i];
								denom_b += v2[i] * v2[i];
							}

							if (denom_a == 0.0 || denom_b == 0.0) {
								return 1.0;
							}

							return dot / (sqrt(denom_a) * sqrt(denom_b));
						};

						double displacement_factor = (
							abs(pt[0]) * std::log(1 + abs(pt[0]))
							+ abs(pt[1]) * std::log(1 + abs(pt[1]))
							) / (2 * MAX_DISP);

						double same_direction_factor = 1.0 - abs(cosine_similarity(pt, major_vector));
						double fitness_factor = cached_objective({ pt[0], pt[1] });

						double a = 0.3, b = 0.8, c = 0.01;
						return (a * fitness_factor* fitness_factor + b * same_direction_factor * displacement_factor + c * displacement_factor) / (a + b + c);
					};

					auto ND_adpated_optimizer = [&MAX_DISP, &previous_step_disp_xy, &major_vector](std::function<double(ND::PT)>& f, ND::PT left, ND::PT right, double tol, int iter)->ND::PT {
						optim::algo_settings_t c_settings;
						c_settings.iter_max = iter;
						//c_settings.rel_sol_change_tol = 0.1;
						c_settings.vals_bound = true; 
						c_settings.lower_bounds = Eigen::VectorXd(2); c_settings.lower_bounds << left[0], left[1];
						c_settings.upper_bounds = Eigen::VectorXd(2); c_settings.upper_bounds << right[0], right[1];

						optim::Mat_t simplex_points(3, 2);
						simplex_points.row(0) << c_settings.lower_bounds(0), c_settings.lower_bounds(1);
						simplex_points.row(1) << c_settings.lower_bounds(0), c_settings.upper_bounds(1);
						simplex_points.row(2) << c_settings.upper_bounds(0), c_settings.upper_bounds(1);
						c_settings.nm_settings.custom_initial_simplex = true; c_settings.nm_settings.initial_simplex_points = simplex_points;
						//c_settings.print_level = 3;

						Eigen::VectorXd c_optimal_displacement(2); c_optimal_displacement << 0.0, 0.0;
						auto ND_adapted_objective_fn = [&f, &MAX_DISP, &previous_step_disp_xy, &major_vector](const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)->double {
							auto cosine_similarity = [](const ND::PT& v1, const ND::PT& v2) {
								double dot = 0.0, denom_a = 0.0, denom_b = 0.0;

								for (unsigned int i = 0; i < v1.size(); ++i) {
									dot += v1[i] * v2[i];
									denom_a += v1[i] * v1[i];
									denom_b += v2[i] * v2[i];
								}

								if (denom_a == 0.0 || denom_b == 0.0) {
									return 1.0;
								}

								return dot / (sqrt(denom_a) * sqrt(denom_b));
							};

							double displacement_factor = (
								abs(vals_inp(0)) * std::log(1 + abs(vals_inp(0)))
								+ abs(vals_inp(1)) * std::log(1 + abs(vals_inp(1)))
								) / (2 * MAX_DISP);

							double same_direction_factor = 1.0 - abs(cosine_similarity({ vals_inp(0), vals_inp(1) }, major_vector));
							double fitness_factor = f({ vals_inp(0), vals_inp(1) });

							double a = 0.3, b = 0.8, c = 0.3;
							return (a * fitness_factor + b * same_direction_factor + c * displacement_factor) / (a + b + c);
						};
						bool success = optim::nm(c_optimal_displacement, ND_adapted_objective_fn, nullptr, c_settings);
						if (true)
							return { c_optimal_displacement(0), c_optimal_displacement(1) };
						else {
							std::cout << "failed suboptimization!" << std::endl;
							return { 0,0 };
						}
					};
					
					auto PSO_adapted_optimizer = [&MAX_DISP, &previous_step_disp_xy](std::function<double(ND::PT)>& f, ND::PT left, ND::PT right, double tol, int iter) {
						return ND::psoOptimization(f, 10, 2, left, right, 10);
					};

					optim::algo_settings_t c_settings;
					c_settings.iter_max = 200;
					//c_settings.rel_objfn_change_tol = 5;
					c_settings.vals_bound = true; c_settings.lower_bounds = -MAX_DISP * Eigen::VectorXd::Ones(2) + previous_step_disp_xy; c_settings.upper_bounds = MAX_DISP * Eigen::VectorXd::Ones(2) + previous_step_disp_xy;

					auto optimal_displacement_vec = ND::custom_splitting_search_nd(ND_adaptetd_f, PSO_adapted_optimizer,
						{ c_settings.lower_bounds(0), c_settings.lower_bounds(1) },
						{ c_settings.upper_bounds(0), c_settings.upper_bounds(1) },
						0.01, 200, 10
					);
					optimal_displacement << optimal_displacement_vec[0], optimal_displacement_vec[1];

					double optimal_fitness = cached_objective({ optimal_displacement(0), optimal_displacement(1) });
					double previous_optimal_fitness = cached_objective({ previous_step_disp_xy(0), previous_step_disp_xy(1) });

					//update major_vector
					if(true)//major_vector[0]==0.0 && major_vector[1]==0.0)
					{
						double vec_x = (1 - optimal_fitness) * respective_polygons.size() * optimal_displacement_vec[0];
						double vec_y = (1 - optimal_fitness) * respective_polygons.size() * optimal_displacement_vec[1];
						double magnitude = sqrt(vec_x * vec_x + vec_y * vec_y);
						if (magnitude != 0 && magnitude== magnitude) {
							major_vector[0]+= vec_x / magnitude; major_vector[1] += vec_y / magnitude;
						}
					}
					std::cout << "major_vector: " << major_vector[0] << ", " << major_vector[1] << std::endl;

					// ignore single geometry alignment if optimal_fitness is over a threshhold
					double single_geom_fitt_threshhold = 0.1;
					if (components_polygons_map[comp_idx].size() == 1 && optimal_fitness > single_geom_fitt_threshhold) {
						optimal_displacement = previous_step_disp_xy;
						optimal_fitness = previous_optimal_fitness;
					}

					// Ignore updating if fitness is high
					if (optimal_fitness > previous_optimal_fitness) {
						optimal_displacement = previous_step_disp_xy;
						optimal_fitness = previous_optimal_fitness;
					}

					

					// adding current component geometries to the temporary components geovector
					for (auto c_poly_ptr : respective_polygons) {

						double ddx = optimal_displacement(0);
						double ddy = optimal_displacement(1);
						bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(ddx, ddy);
						auto translated_geometry = translate_geometry(c_poly_ptr->get_definition(), trans_obj);
						auto stitched_pixels = RPR.readPolygonPixels<float>(translated_geometry, RasterPixelsStitcherStartegy::contours);
						auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
						double c_polygon_fitness = fitness_from_stats_functor(stats);
						double c_polygon_previous_fitness = c_poly_ptr->get_double_attribute(FITT_COLUMN_NAME);
						double c_polygon_previous_disparity_x = c_poly_ptr->get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first);
						double c_polygon_previous_disparity_y = c_poly_ptr->get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second);

						bool alignment_success = c_polygon_fitness < c_polygon_previous_fitness;

						// add to component temporary shp
						if (!alignment_success) {
							bg::strategy::transform::translate_transformer<double, 2, 2> prev_trans_obj(
								c_polygon_previous_disparity_x,
								c_polygon_previous_disparity_y
							);
							component_gvector.add_geometry(translate_geometry(c_poly_ptr->get_definition(), prev_trans_obj));
						}
						else
							component_gvector.add_geometry(translated_geometry);

						auto saved_gwa = component_gvector.geometries_container.rbegin();
						saved_gwa->set_double_attribute("obj", optimal_fitness);
						saved_gwa->set_int_attribute("cmp_id", comp_idx);
						saved_gwa->set_double_attribute("dispx", (alignment_success) ? ddx : c_polygon_previous_disparity_x);
						saved_gwa->set_double_attribute("dispy", (alignment_success) ? ddy : c_polygon_previous_disparity_y);
						saved_gwa->set_double_attribute("ind_obj", (alignment_success) ? c_polygon_fitness : c_polygon_previous_fitness);

						// updating fields
						c_poly_ptr->set_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first, optimal_displacement(0));
						c_poly_ptr->set_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second, optimal_displacement(1));
					}

				}
				bar.finish();
				std::string comp_file = params->temp_dir + "//comp_" + std::to_string(*distance_val_iter) + ".shp";
				component_gvector.to_file(comp_file);
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
			auto linestring_fitness_evaluator = [&RPR, &null_value, &fitness_from_stats_functor](const Boost_LineString_2& ref_geometry)->double {
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
					

					std::function<double(double)> bound_objective = [&RPR, &r2r_constants, &null_value, &linestring_fitness_evaluator, &component_exterior_edges](double c_component_height)->double {
						double total_fitness = 0;
						double total_length = 0;
						for (auto& c_exterior_edge : component_exterior_edges) {
							double c_edge_length = bg::length(c_exterior_edge);
							bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(r2r_constants.first * c_component_height, r2r_constants.second * c_component_height);
							auto translated_geometry = translate_geometry(c_exterior_edge, trans_obj);
							double c_geometry_obj_val = linestring_fitness_evaluator(translated_geometry);
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
					double previous_optimal_fitness = cached_objective(0.0); // TODO:: check this (cached objective applies height on original geometries) maybe should be previous_step_disp_value

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