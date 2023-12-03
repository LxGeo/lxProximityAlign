#include "nm_search.h"
#include "boost/multi_array.hpp"

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
					
					major_vector = { previous_step_disp_x_value, previous_step_disp_y_value };
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

						double a = 0.5, b = 0.5, c = 0.01;
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
						return ND::psoOptimization(f, 5, 2, left, right, 10);
					};

					optim::algo_settings_t c_settings;
					c_settings.iter_max = 100;
					//c_settings.rel_objfn_change_tol = 5;
					c_settings.vals_bound = true; c_settings.lower_bounds = -MAX_DISP * Eigen::VectorXd::Ones(2) + previous_step_disp_xy; c_settings.upper_bounds = MAX_DISP * Eigen::VectorXd::Ones(2) + previous_step_disp_xy;

					auto optimal_displacement_vec = ND::custom_splitting_search_nd(ND_adaptetd_f, PSO_adapted_optimizer,
						{ c_settings.lower_bounds(0), c_settings.lower_bounds(1) },
						{ c_settings.upper_bounds(0), c_settings.upper_bounds(1) },
						0.01, c_settings.iter_max, 4
					);
					optimal_displacement << optimal_displacement_vec[0], optimal_displacement_vec[1];

					double optimal_fitness = cached_objective({ optimal_displacement(0), optimal_displacement(1) });
					double previous_optimal_fitness = cached_objective({ previous_step_disp_xy(0), previous_step_disp_xy(1) });

					//update major_vector
					//update only in first iteration
					/*
					if (distance_val_iter == neighbour_distance_band_values.rbegin()) {
						double magnitude = sqrt(optimal_displacement_vec[0] * optimal_displacement_vec[0] + optimal_displacement_vec[1] * optimal_displacement_vec[1]);
						if (magnitude != 0 && magnitude == magnitude) {
							major_vector[0] += (1 - optimal_fitness) * optimal_displacement_vec[0] / magnitude;
							major_vector[1] += (1 - optimal_fitness) * optimal_displacement_vec[1] / magnitude;
						}
					}
					*/

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
						double c_polygon_previous_disparity_x = c_poly_ptr->get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first);
						double c_polygon_previous_disparity_y = c_poly_ptr->get_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second);
						double c_polygon_previous_fitness = c_poly_ptr->get_double_attribute(FITT_COLUMN_NAME);

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
						c_poly_ptr->set_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.first, (alignment_success) ? ddx : c_polygon_previous_disparity_x);
						c_poly_ptr->set_double_attribute(OBJECTIVE_FIELD_NAME_PAIR.second, (alignment_success) ? ddy : c_polygon_previous_disparity_y);
						c_poly_ptr->set_double_attribute(FITT_COLUMN_NAME, (alignment_success) ? c_polygon_fitness : c_polygon_previous_fitness);
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
					//if (c_fitness_value < FITNESS_THRESH)
					//	continue;

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
					std::list<Boost_LineString_2> component_exterior_edges;

					/*
					// Extract component exterior edges using arrangment
					Arrangement respective_polygons_arrangment = ArrangmentFromPolygons(respective_polygons);
										
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
					*/
					for (auto& c_poly : respective_polygons) {
						auto& c_geom = c_poly->get_definition();
						for (int edge_idx = 0; edge_idx < c_geom.outer().size() - 1; edge_idx++) {
							component_exterior_edges.push_back(
								Boost_LineString_2({ c_geom.outer().at(edge_idx), c_geom.outer().at(edge_idx+1)  })
							);
						}
					}

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


		void nm_proximity_align_linear(
			std::unordered_map<std::string, matrix>& matrices_map, GeoImage<cv::Mat>& ref_gimg,
			GeoVector<Boost_LineString_2>& input_geovector,
			std::vector<double>& neighbour_distance_band_values,
			std::function<float(numcpp::DetailedStats<float>&)> fitness_from_stats_functor,
			double MAX_DISP,
			std::pair<std::string, std::string> OBJECTIVE_FIELD_NAME_PAIR) {


			SpatialWeights<Boost_LineString_2> PSW = SpatialWeights<Boost_LineString_2>::from_geovector(input_geovector);
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

			auto linestring2seg_unpooling_fn = [](const Geometries_with_attributes<Boost_LineString_2>& parent_linestring) -> std::list<Geometries_with_attributes<Boost_LineString_2>> {
				std::list<Geometries_with_attributes<Boost_LineString_2>> unpooled_linestrings;

				auto& c_linestring = parent_linestring.get_definition();
				auto c_line_vertex_it = c_linestring.begin();
				auto prev = c_line_vertex_it;
				int child_position = 0;
				for (++c_line_vertex_it; c_line_vertex_it != c_linestring.end(); ++c_line_vertex_it, ++child_position) {
					Boost_LineString_2 c_segment;
					c_segment.push_back(*prev);
					c_segment.push_back(*c_line_vertex_it);
					Geometries_with_attributes<Boost_LineString_2> c_gwa(c_segment);
					c_gwa.set_int_attribute("pid", parent_linestring.get_int_attribute(GeoVector<Boost_LineString_2>::ID_FIELD_NAME));
					c_gwa.set_int_attribute("pos", child_position);
					unpooled_linestrings.push_back(c_gwa);
					prev = c_line_vertex_it;
				}
				return unpooled_linestrings;
			};

			std::function<Geometries_with_attributes<Boost_LineString_2>(std::list<Geometries_with_attributes<Boost_LineString_2>>&)> reconstruction_fn = [](std::list<Geometries_with_attributes<Boost_LineString_2>>& geom_parts) ->Geometries_with_attributes<Boost_LineString_2> {
				geom_parts.sort([](const Geometries_with_attributes<Boost_LineString_2>& a, const Geometries_with_attributes<Boost_LineString_2>& b) { return a.get_double_attribute("pos") < b.get_double_attribute("pos"); });
				Boost_LineString_2 reconstructed_l;
				if (geom_parts.size() == 0) { return Geometries_with_attributes<Boost_LineString_2>(); }
				if (geom_parts.size() == 1) {
					return Geometries_with_attributes<Boost_LineString_2>(geom_parts.begin()->get_definition());
				}
				auto first_part = geom_parts.begin()->get_definition();
				auto second_part = std::next(geom_parts.begin())->get_definition();
				if (bg::distance(first_part.at(0), second_part.at(0)) < 1e-5 || bg::distance(first_part.at(0), second_part.at(1)) < 1e-5)
					reconstructed_l.push_back(geom_parts.begin()->get_definition().at(1)); // adding starting point
				else
					reconstructed_l.push_back(geom_parts.begin()->get_definition().at(0)); // adding starting point
				for (const auto& geom_part : geom_parts) {
					auto point_to_add = (bg::distance(reconstructed_l.back(), geom_part.get_definition().at(0)) < 1e-5) ? geom_part.get_definition().at(1) : geom_part.get_definition().at(0);
					reconstructed_l.push_back(point_to_add);
				};
				Geometries_with_attributes<Boost_LineString_2> gwa_l(reconstructed_l);
				return gwa_l;
			};

			GeoVector<Boost_LineString_2> unpooled_linestring_gvec; unpool_geovector<Boost_LineString_2, Boost_LineString_2>(input_geovector, unpooled_linestring_gvec, linestring2seg_unpooling_fn);

			std::function<LinearTopology<EK>::SegmentIdentification(const Geometries_with_attributes<Boost_LineString_2>&)> gwa2segId_fn =
				[](const Geometries_with_attributes<Boost_LineString_2>& gwa) ->LinearTopology<EK>::SegmentIdentification {
				LinearTopology<EK>::SegmentIdentification seg_id; seg_id.parent_id = gwa.get_int_attribute("pid"); seg_id.position = gwa.get_int_attribute("pos");
				return seg_id;
			};
			LinearTopology<EK>::Arrangement_2 arr = LinearTopology<EK>::arrangmentFromLineStringGeovector(unpooled_linestring_gvec, gwa2segId_fn);
			// Defining datastructures
			std::unordered_map<LinearTopology<EK>::Arrangement_2::Vertex_const_iterator, std::pair<double, double>> displacement_map;
			for (auto v_handle = arr.vertices_begin(); v_handle != arr.vertices_end(); v_handle++) {
				displacement_map[v_handle] = { 0,0 };
			}

			/*for (int c_iter = 0; c_iter < 10; c_iter++) {

				std::pair<double, double> c_global_disp = { 0.0, 0.0 };
				// Loop align each vertex indep
				GeoVector<Boost_Point_2> pts_gvec;
				for (auto c_vertex = arr.vertices_begin(); c_vertex != arr.vertices_end(); c_vertex++) {

					std::function<double(const std::pair<double, double>&)> bound_objective = [&](const std::pair<double, double>& disp) ->double {

						// get all respective segments after last modification
						std::list<Boost_LineString_2> respective_segments;
						LinearTopology<EK>::Arrangement_2::Halfedge_around_vertex_const_circulator  circ_first, circ_current;
						circ_first = circ_current = c_vertex->incident_halfedges();
						do {
							Boost_LineString_2 c_incident(
								{
									{PointTraits<Point_2>::getX(c_vertex->point()) + displacement_map[c_vertex].first, PointTraits<Point_2>::getY(c_vertex->point()) + displacement_map[c_vertex].second},
									{PointTraits<Point_2>::getX(circ_current->source()->point()) + displacement_map[circ_current->source()].first, PointTraits<Point_2>::getY(circ_current->source()->point()) + displacement_map[circ_current->source()].second}
								}
							);
							respective_segments.push_back(c_incident);
						} while (++circ_current != circ_first);

						// transform edges
						std::list<Boost_LineString_2> transformed_respective_segments;
						for (auto& c_geom : respective_segments) {
							Boost_LineString_2 transformed_geom(c_geom);
							auto& src_pt = c_geom.at(0);
							auto& tar_pt = c_geom.at(1);
							transformed_geom.at(0).set<0>(src_pt.get<0>() + disp.first); transformed_geom.at(0).set<1>(src_pt.get<1>() + disp.second);
							transformed_geom.at(1).set<0>(tar_pt.get<0>() + disp.first); transformed_geom.at(1).set<1>(tar_pt.get<1>() + disp.second);
							transformed_respective_segments.push_back(transformed_geom);
						}

						// compute fitness
						double total_fitness = 0.0, total_length = 0.0;
						for (auto& c_geom : transformed_respective_segments) {
							double c_geometry_obj_val = linestring_fitness_evaluator(c_geom);
							double c_edge_length = bg::length(c_geom);
							total_fitness += (c_geometry_obj_val * c_edge_length);
							total_length += c_edge_length;
						}
						double mean_fitness = total_fitness / total_length;
						return mean_fitness;
					};

					auto ND_adpated_optimizer = [&MAX_DISP](std::function<double(const std::pair<double, double>&)>& f, ND::PT left, ND::PT right, double tol, int iter)->ND::PT {
						optim::algo_settings_t c_settings;
						c_settings.iter_max = iter;
						//c_settings.rel_sol_change_tol = 0.1;
						c_settings.vals_bound = true;
						c_settings.lower_bounds = Eigen::VectorXd(2); c_settings.lower_bounds << left[0], left[1];
						c_settings.upper_bounds = Eigen::VectorXd(2); c_settings.upper_bounds << right[0], right[1];

						optim::Mat_t simplex_points(3, 2);
						simplex_points.row(0) << -5, 10;
						simplex_points.row(1) << 0, 15;
						simplex_points.row(2) << -10, -5;
						c_settings.nm_settings.custom_initial_simplex = true; c_settings.nm_settings.initial_simplex_points = simplex_points;
						//c_settings.print_level = 3;

						Eigen::VectorXd c_optimal_displacement(2); c_optimal_displacement << 0.0, 0.0;
						auto ND_adapted_objective_fn = [&f, &MAX_DISP](const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)->double {
							double fitness_factor = f({ vals_inp(0), vals_inp(1) });
							return fitness_factor;
						};
						bool success = optim::nm(c_optimal_displacement, ND_adapted_objective_fn, nullptr, c_settings);
						if (true)
							return { c_optimal_displacement(0), c_optimal_displacement(1) };
						else {
							std::cout << "failed suboptimization!" << std::endl;
							return { 0,0 };
						}
					};

					auto optimal_displacment = ND_adpated_optimizer(bound_objective, { -20,-20 }, { 20,20 }, 0, 50);
					//c_global_disp.first += optimal_displacment[0]/ arr.number_of_vertices(); c_global_disp.second += optimal_displacment[1]/ arr.number_of_vertices();
					displacement_map[c_vertex] = { optimal_displacment[0], optimal_displacment[1] };
					Geometries_with_attributes<Boost_Point_2> c_pt(Boost_Point_2(PointTraits<Point_2>::getX(c_vertex->point()), PointTraits<Point_2>::getY(c_vertex->point())));
					c_pt.set_double_attribute("dx", optimal_displacment[0]);
					c_pt.set_double_attribute("dy", optimal_displacment[1]);
					c_pt.set_double_attribute("fitt", bound_objective({ optimal_displacment[0], optimal_displacment[1]}));
					pts_gvec.add_geometry(c_pt);
				}
				pts_gvec.to_file(params->temp_dir + "/vertices_disp.shp");

				for (auto c_vertex = arr.vertices_begin(); c_vertex != arr.vertices_end(); c_vertex++) {
					//displacement_map[c_vertex] = c_global_disp;
				}
				

				auto aligned_iter_gvec = LinearTopology<EK>::FilteredDisplacedGeovectorFromArrangment(arr, [&](LinearTopology<EK>::Arrangement_2::Vertex_const_iterator viter) {return true; }, displacement_map);
				GeoVector<Boost_LineString_2> pooled_gvec;
				pool_geovector(aligned_iter_gvec, pooled_gvec, reconstruction_fn);
				pooled_gvec.to_file(params->temp_dir + "/_"+std::to_string(c_iter) + "_aligned_full.shp");

			}*/

			auto ND_adpated_optimizer = [&MAX_DISP](std::function<double(const std::pair<double, double>&)>& f, ND::PT left, ND::PT right, double tol, int iter)->ND::PT {
				optim::algo_settings_t c_settings;
				c_settings.iter_max = iter;
				//c_settings.rel_sol_change_tol = 0.1;
				c_settings.vals_bound = true;
				c_settings.lower_bounds = Eigen::VectorXd(2); c_settings.lower_bounds << left[0], left[1];
				c_settings.upper_bounds = Eigen::VectorXd(2); c_settings.upper_bounds << right[0], right[1];

				optim::Mat_t simplex_points(3, 2);
				/*simplex_points.row(0) << -5, 10;
				simplex_points.row(1) << 0, 15;
				simplex_points.row(2) << -10, -5;*/
				simplex_points.row(0) << -50, 50;
				simplex_points.row(1) << 50, 50;
				simplex_points.row(2) << -50, -50;
				c_settings.nm_settings.custom_initial_simplex = true; c_settings.nm_settings.initial_simplex_points = simplex_points;
				//c_settings.print_level = 3;

				Eigen::VectorXd c_optimal_displacement(2); c_optimal_displacement << 0.0, 0.0;
				auto ND_adapted_objective_fn = [&f, &MAX_DISP](const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)->double {
					double fitness_factor = f({ vals_inp(0), vals_inp(1) });
					return fitness_factor;
				};
				bool success = optim::nm(c_optimal_displacement, ND_adapted_objective_fn, nullptr, c_settings);
				std::cout << c_optimal_displacement(0) << " " << c_optimal_displacement(1) << std::endl;
				if (true)
					return { c_optimal_displacement(0), c_optimal_displacement(1) };
				else {
					std::cout << "failed suboptimization!" << std::endl;
					return { 0,0 };
				}
			};

			for (auto distance_val_iter = neighbour_distance_band_values.rbegin(); distance_val_iter != neighbour_distance_band_values.rend(); ++distance_val_iter) {
				auto disconnection_lambda = [&distance_val_iter](double x)->bool {return x > *distance_val_iter; };
				PSW.disconnect_edges(disconnection_lambda);
				PSW.run_labeling();
				std::cout << "Aligning N components: " << PSW.n_components << std::endl;

				std::map<size_t, std::unordered_set<size_t>> components_geometries_map;
				std::map<size_t, Eigen::VectorXd> components_transforms_map;
				for (size_t _idx = 0; _idx < PSW.component_labels.size(); ++_idx) {
					size_t comp_id = PSW.component_labels[_idx];
					int geom_id = PSW.geometries_container[_idx].get_int_attribute(GeoVector<Boost_LineString_2>::ID_FIELD_NAME);
					if (components_geometries_map.find(comp_id) != components_geometries_map.end())
						components_geometries_map[comp_id].insert(geom_id);
					else
					{
						components_geometries_map[comp_id] = std::unordered_set<size_t>();
						components_geometries_map[comp_id].insert(geom_id);
					}
				}

				std::unordered_map<size_t, std::set<LinearTopology<EK>::Arrangement_2::Vertex_const_iterator>> components_verices_map;
				// Fill cluster vertices map with verices correpsonding to each component
				for (auto v_handle = arr.vertices_begin(); v_handle != arr.vertices_end(); v_handle++) {
					LinearTopology<EK>::Arrangement_2::Halfedge_around_vertex_const_circulator  circ_first, circ_current;
					circ_first = circ_current = v_handle->incident_halfedges();
					do {
						auto pid = circ_current->curve().data().front().parent_id;
						for (auto& [component_id, respective_geometries_set] : components_geometries_map) {
							if (respective_geometries_set.find(pid) != respective_geometries_set.end()) {
								components_verices_map[component_id].insert(v_handle);
							}
						}
					}
					while (++circ_current != circ_first);
				}

				GeoVector<Boost_Point_2> pts_gvec;
				for (const auto& [cluster_id, respective_vertices_handles] : components_verices_map) {

					std::function<bool(LinearTopology<EK>::Arrangement_2::Vertex_const_iterator)> filter_predicate = [&](LinearTopology<EK>::Arrangement_2::Vertex_const_iterator viter) {
						return respective_vertices_handles.find(viter) != respective_vertices_handles.end();
					};

					std::function<double(const std::pair<double, double>&)> bound_objective = [&](const std::pair<double, double>& disp) ->double {						
						// TODO apply only on respective 
						std::unordered_map < LinearTopology<EK>::Arrangement_2::Vertex_const_iterator, std::pair<double, double> > c_cluster_disp_map;
						for (auto& c_v_handle : respective_vertices_handles)
							c_cluster_disp_map[c_v_handle] = { displacement_map[c_v_handle].first + disp.first, displacement_map[c_v_handle].second + disp.second };
						auto filtered_gvec = LinearTopology<EK>::FilteredDisplacedGeovectorFromArrangment(arr, filter_predicate, c_cluster_disp_map);

						double total_fitness = 0.0, total_length = 0.0;
						for (auto& c_gwa : filtered_gvec) {
							double c_geometry_obj_val = linestring_fitness_evaluator(c_gwa.get_definition());
							double c_edge_length = bg::length(c_gwa.get_definition());
							total_fitness += (c_geometry_obj_val * c_edge_length);
							total_length += c_edge_length;
						}
						double mean_fitness = total_fitness / total_length;
						return mean_fitness;
					};

					auto optimal_displacment = ND_adpated_optimizer(bound_objective, { -MAX_DISP,-MAX_DISP }, { MAX_DISP,MAX_DISP }, 0, 100);

					for (auto& c_v_handle : respective_vertices_handles) {
						displacement_map[c_v_handle].first += optimal_displacment[0];
						displacement_map[c_v_handle].second += optimal_displacment[1];
					}

					for (auto v_handle : respective_vertices_handles) {
						Geometries_with_attributes<Boost_Point_2> c_pt(Boost_Point_2(PointTraits<Point_2>::getX(v_handle->point()), PointTraits<Point_2>::getY(v_handle->point())));
						c_pt.set_double_attribute("dx", optimal_displacment[0]);
						c_pt.set_double_attribute("dy", optimal_displacment[1]);
						pts_gvec.add_geometry(c_pt);
					}
				}
				pts_gvec.to_file(params->temp_dir + "/_" + std::to_string(*distance_val_iter) + "vertices_disp.shp");				
			}
			
			// per vertex refinement
			/*for (auto v_handle = arr.vertices_begin(); v_handle != arr.vertices_end(); v_handle++) {

				std::function<double(const std::pair<double, double>&)> bound_objective = [&](const std::pair<double, double>& disp) ->double {
					// TODO apply only on respective 
					std::unordered_map<LinearTopology<EK>::Arrangement_2::Vertex_const_iterator, std::pair<double, double> > c_cluster_disp_map;
					c_cluster_disp_map[v_handle] = { displacement_map[v_handle].first + disp.first, displacement_map[v_handle].second + disp.second };
					LinearTopology<EK>::Arrangement_2::Halfedge_around_vertex_const_circulator  circ_first, circ_current;
					circ_first = circ_current = v_handle->incident_halfedges();
					do {
						c_cluster_disp_map[circ_current->source()] = { 0,0 };
					} while (++circ_current != circ_first);
					
					auto filtered_gvec = LinearTopology<EK>::FilteredDisplacedGeovectorFromArrangment(arr, [&v_handle](LinearTopology<EK>::Arrangement_2::Vertex_const_iterator v) {return v == v_handle; }, c_cluster_disp_map);

					double total_fitness = 0.0, total_length = 0.0;
					for (auto& c_gwa : filtered_gvec) {
						double c_geometry_obj_val = linestring_fitness_evaluator(c_gwa.get_definition());
						double c_edge_length = bg::length(c_gwa.get_definition());
						total_fitness += (c_geometry_obj_val * c_edge_length);
						total_length += c_edge_length;
					}
					double mean_fitness = total_fitness / total_length;
					return mean_fitness;
				};
				auto optimal_displacment = ND_adpated_optimizer(bound_objective, { -5,-5 }, { 5,5 }, 0, 50);
				displacement_map[v_handle].first += optimal_displacment[0];
				displacement_map[v_handle].second += optimal_displacment[1];

			}*/

			auto aligned_gvec = LinearTopology<EK>::FilteredDisplacedGeovectorFromArrangment(arr, [&](LinearTopology<EK>::Arrangement_2::Vertex_const_iterator viter) {return true; }, displacement_map);
			//aligned_gvec.to_file("C:/DATA_SANDBOX/fault_to_align/data_mo_TIZ110723/data_mo_TIZ110723/example1Bishop/output_dir/aligned_parts.shp");

			GeoVector<Boost_LineString_2> pooled_gvec;
			pool_geovector(aligned_gvec, pooled_gvec, reconstruction_fn);
			pooled_gvec.to_file(params->temp_dir + "/aligned_full.shp");

 		}

	}
}