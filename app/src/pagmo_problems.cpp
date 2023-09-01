#include "pagmo_problems.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{

		template <typename range>
		Eigen::MatrixXd computePairwiseDistance(const range& points) {
			int N = points.size();

			// Convert the set of points to an Eigen matrix
			Eigen::MatrixXd data(N, 2);
			size_t idx = 0;
			for (const auto& p : points) {
				data(idx, 0) = PointTraits<std::ranges::range_value_t<range>>::getX(p);
				data(idx, 1) = PointTraits<std::ranges::range_value_t<range>>::getY(p);
				idx++;
			}

			// Compute squared norms for each row in data
			Eigen::VectorXd squaredNorms = data.rowwise().squaredNorm();

			// Use broadcasting to compute the pairwise squared distances
			Eigen::MatrixXd pairwiseSquaredDistances = squaredNorms.rowwise().replicate(data.rows())
				+ squaredNorms.transpose().colwise().replicate(data.rows())
				- 2 * data * data.transpose();

			// Ensure no negative values due to numerical inaccuracies
			pairwiseSquaredDistances = pairwiseSquaredDistances.cwiseMax(0.0);

			// Compute the actual distances
			Eigen::MatrixXd distances = pairwiseSquaredDistances.array().sqrt();
			return distances;
		}

		template <typename range>
		Eigen::MatrixXd computePairwiseDifference(const range& points) {
			int N = points.size();

			// Convert the set of points to an Eigen matrix
			Eigen::MatrixXd data(N, 2);
			size_t idx = 0;
			for (const auto& p : points) {
				data(idx, 0) = PointTraits<std::ranges::range_value_t<range>>::getX(p);
				data(idx, 1) = PointTraits<std::ranges::range_value_t<range>>::getY(p);
				idx++;
			}

			// Create a 3D tensor to store pairwise differences
			// For simplicity, we'll use a MatrixXd and reshape it later
			Eigen::MatrixXd pairwiseDifferences(N * N, 2);

			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					pairwiseDifferences.row(i * N + j) = data.row(i) - data.row(j);
				}
			}

			return pairwiseDifferences;
		}

		Point_2 rotate_around_point(const Point_2& A, const Point_2& B, double angle_in_degrees) {
			// Convert angle to radians
			double angle_in_radians = angle_in_degrees * CGAL_PI / 180.0;

			auto cos_angle = std::cos(angle_in_radians);
			auto sin_angle = std::sin(angle_in_radians);

			// Translate B so that A becomes the origin
			auto x_translated = B.x() - A.x();
			auto y_translated = B.y() - A.y();

			// Rotate around the origin
			auto x_rotated = x_translated * cos_angle - y_translated * sin_angle;
			auto y_rotated = x_translated * sin_angle + y_translated * cos_angle;

			// Translate back
			auto x_final = x_rotated + A.x();
			auto y_final = y_rotated + A.y();

			return Point_2(x_final, y_final);
		}

		void pagmo_proximity_align_linear(
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
			std::unordered_map<LinearTopology<EK>::Arrangement_2::Vertex_const_iterator, std::list<std::tuple<double, double, double>>> displacement_map;
			
			tqdm bar;
			int c_estimated_vertex_handle = 0;
			for (auto v_handle = arr.vertices_begin(); v_handle != arr.vertices_end(); v_handle++) {
				bar.progress(c_estimated_vertex_handle++, arr.number_of_vertices());
				size_t N_MAX_LEVEL = 1;
				auto temp_half_edge = v_handle->incident_halfedges()->source()->incident_halfedges();
				while (temp_half_edge->source() != v_handle) 
					++temp_half_edge;
				std::map<size_t, std::set<LinearTopology<EK>::Arrangement_2::Halfedge_around_vertex_const_circulator>> N_neighbours_edges;
				N_neighbours_edges[0].insert(temp_half_edge);
				std::set<LinearTopology<EK>::Arrangement_2::Vertex_const_iterator> all_neighbour_vertices = {v_handle,};
				std::map<size_t, std::set<LinearTopology<EK>::Arrangement_2::Vertex_const_iterator>> N_neighbours_vertices = { {0, {v_handle,}} };
				// filling datastructures
				for (int c_neighbour_level = 1; c_neighbour_level <= N_MAX_LEVEL; c_neighbour_level++) {
					//for previous level elements
					for (auto c_prev_neighbour_edge : N_neighbours_edges[c_neighbour_level - 1]) {

						LinearTopology<EK>::Arrangement_2::Halfedge_around_vertex_const_circulator  circ_first, circ_current;
						circ_first = circ_current = c_prev_neighbour_edge->source()->incident_halfedges();
						do {
							bool s_already_added = all_neighbour_vertices.find(circ_current->source()) != all_neighbour_vertices.end();
							bool t_already_added = all_neighbour_vertices.find(circ_current->target()) != all_neighbour_vertices.end();
							bool already_added = s_already_added && t_already_added;
							
							if (!already_added) {
								N_neighbours_edges[c_neighbour_level].insert(circ_current);
								N_neighbours_vertices[c_neighbour_level].insert(circ_current->source());
								all_neighbour_vertices.insert(circ_current->source());
								all_neighbour_vertices.insert(circ_current->target());
							}
						} while (++circ_current != circ_first);
					}
				}
				N_neighbours_edges[0].clear();

				
				/*cluster_problem<int>::bound_type objective_boundary = {
						std::vector<double>(2* all_neighbour_vertices.size(), -MAX_DISP),
						std::vector<double>(2* all_neighbour_vertices.size(), MAX_DISP)
				};*/
				cluster_problem<int>::bound_type objective_boundary;
				double MAX_ROT_DEGREES = 10.0;
				objective_boundary.first = std::vector<double>(all_neighbour_vertices.size() + 1 ,-MAX_ROT_DEGREES);
				objective_boundary.second = std::vector<double>(all_neighbour_vertices.size() + 1, MAX_ROT_DEGREES);

				
				auto all_neighbour_points = all_neighbour_vertices | std::views::transform([](auto& v_handle) { return v_handle->point(); });
				std::function<double(const std::vector<double>&)> bound_objective = [&](const std::vector<double>& disp) ->double {					
					double tx = disp[0] / MAX_ROT_DEGREES * MAX_DISP, ty = disp[1] / MAX_ROT_DEGREES * MAX_DISP;
					// disp contain displacement values for each vertex in the following order [ v0_x, v0_y, v1_x, ..., vN_x, vN_y]
					std::unordered_map < LinearTopology<EK>::Arrangement_2::Vertex_const_iterator, std::pair<double, double> > c_neighbours_disp_map;
					size_t c_vertex_idx = 0;
					Point_2 displaced_center = v_handle->point() + Vector_2(tx, ty);
					for (auto& c_v_handle : all_neighbour_vertices) {
						if (c_v_handle == v_handle) {
							c_neighbours_disp_map[c_v_handle] = { tx,  ty };
						}
						else {
							Point_2 new_b_position = rotate_around_point(displaced_center, c_v_handle->point() + Vector_2(tx, ty), disp[2 + c_vertex_idx]);
							auto diff_combined_transfrom = new_b_position - c_v_handle->point();							
							c_neighbours_disp_map[c_v_handle] = { CGAL::to_double(diff_combined_transfrom.x()), CGAL::to_double(diff_combined_transfrom.y()) };
							c_vertex_idx++;
						}						
					}


					std::list<Boost_LineString_2> transformed_edges;
					
					for (const auto& [c_level, c_neighbours_edges_at_level] : N_neighbours_edges) {
						for (auto& c_edge_handle : c_neighbours_edges_at_level) {
							auto source_point = c_edge_handle->source();
							auto target_point = c_edge_handle->target();
							Boost_LineString_2 transformed_edge = {
								{ PointTraits<Point_2>::getX(source_point->point()) + c_neighbours_disp_map[source_point].first, PointTraits<Point_2>::getY(source_point->point()) + c_neighbours_disp_map[source_point].second },
								{ PointTraits<Point_2>::getX(target_point->point()) + c_neighbours_disp_map[target_point].first, PointTraits<Point_2>::getY(target_point->point()) + c_neighbours_disp_map[target_point].second },
							};
							transformed_edges.push_back(transformed_edge);
						}
					}

					double total_fitness = 0.0, total_length = 0.0;
					for (auto& c_gwa : transformed_edges) {
						double c_geometry_obj_val = linestring_fitness_evaluator(c_gwa);
						double c_edge_length = bg::length(c_gwa);
						total_fitness += (c_geometry_obj_val * c_edge_length);
						total_length += c_edge_length;
					}
					double mean_fitness = total_fitness / total_length;
					return mean_fitness;
				};

				pagmo::problem problem{cluster_problem{bound_objective, objective_boundary }};
				pagmo::population pop{ problem, 5 };
				pop.set_x(0, std::vector<double>(objective_boundary.first.size(), 0.0));

				//pagmo::algorithm algo{pagmo::nlopt("neldermead")};
				pagmo::algorithm algo{pagmo::pso{20}};
				pop = algo.evolve(pop);
				std::vector<double> best_value = pop.champion_f();
				std::vector<double> best_arg = pop.champion_x();

				size_t c_vertex_idx = 0;
				double tx = best_arg[0] / MAX_ROT_DEGREES * MAX_DISP, ty = best_arg[1] / MAX_ROT_DEGREES * MAX_DISP;
				Point_2 displaced_center = v_handle->point() + Vector_2(tx, ty);
				for (auto& c_v_handle : all_neighbour_vertices) {
					if (c_v_handle == v_handle) {
						if  (best_value[0]<0.9)
							displacement_map[c_v_handle].push_back({ tx, ty, 1.0 / (0.0001+best_value[0]) });
						//else
						//	displacement_map[c_v_handle].push_back({ 0, 0, 1.0 });
					}
					else {
						if (best_value[0] < 0.9) {
							Point_2 new_b_position = rotate_around_point(displaced_center, c_v_handle->point() + Vector_2(tx, ty), best_arg[2 + c_vertex_idx]);
							auto diff_combined_transfrom = new_b_position - c_v_handle->point();
							//displacement_map[c_v_handle].push_back({ CGAL::to_double(diff_combined_transfrom.x()), CGAL::to_double(diff_combined_transfrom.y()), 1.0 / (0.0001 + best_value[0]) });
						}
						//else
						//	displacement_map[c_v_handle].push_back({ 0,0,1.0 });
						c_vertex_idx++;
					}
				}

			}

			bar.finish();
			std::unordered_map<LinearTopology<EK>::Arrangement_2::Vertex_const_iterator, std::pair<double, double>> agg_displacement_map;
			for (const auto& [k, v] : displacement_map) {
				double x_sum = 0, y_sum = 0, z_sum = 0;
				for (auto& el : v) {
					double x, y, z;
					std::tie(x, y, z) = el;
					x_sum += x*z;
					y_sum += y*z;
					z_sum += z;
				}
				agg_displacement_map[k] = { x_sum / z_sum, y_sum / z_sum };
			}

			/*
			for (auto distance_val_iter = neighbour_distance_band_values.rbegin(); distance_val_iter != neighbour_distance_band_values.rend(); ++distance_val_iter) {
				auto disconnection_lambda = [&distance_val_iter](double x)->bool {return x > *distance_val_iter; };
				PSW.disconnect_edges(disconnection_lambda);
				PSW.run_labeling();
				std::cout << "Aligning N components: " << PSW.n_components << std::endl;

				// Creating cluster based on distance between features
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

				// Creating cluster of vertices from the already clustered geometries
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
					} while (++circ_current != circ_first);
				}

				GeoVector<Boost_Point_2> pts_gvec;
				// Align clusters 1by1
				for (const auto& [cluster_id, respective_vertices_handles] : components_verices_map) {
					std::cout << "Processing cluster: " << cluster_id << " with N vertices: " << respective_vertices_handles.size() << std::endl;
					std::function<bool(LinearTopology<EK>::Arrangement_2::Vertex_const_iterator)> filter_predicate = [&](LinearTopology<EK>::Arrangement_2::Vertex_const_iterator viter) {
						return respective_vertices_handles.find(viter) != respective_vertices_handles.end();
					};

					auto get_unique_vertex_objective = [&](LinearTopology<EK>::Arrangement_2::Vertex_const_iterator unique_vertex) {

						std::function<double(const std::vector<double>&)> bound_objective_all = [&](const std::vector<double>& disp) ->double {
							// TODO apply only on respective 
							std::unordered_map < LinearTopology<EK>::Arrangement_2::Vertex_const_iterator, std::pair<double, double> > c_cluster_disp_map;
							for (auto& c_v_handle : respective_vertices_handles)
							{
								c_cluster_disp_map[c_v_handle] = { displacement_map[c_v_handle].first + disp[2], displacement_map[c_v_handle].second + disp[3] };
							}
							c_cluster_disp_map[unique_vertex].first = disp[0]; c_cluster_disp_map[unique_vertex].second = disp[1];
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
						return bound_objective_all;
					};
					cluster_problem<int>::bound_type objective_boundary_all = { 
						std::vector<double>(4, -MAX_DISP),
						std::vector<double>(4, MAX_DISP) 
					};


					std::function<double(const std::vector<double>&)> bound_objective = [&](const std::vector<double>& disp) ->double {
						// TODO apply only on respective 
						std::unordered_map < LinearTopology<EK>::Arrangement_2::Vertex_const_iterator, std::pair<double, double> > c_cluster_disp_map;
						for (auto& c_v_handle : respective_vertices_handles)
							c_cluster_disp_map[c_v_handle] = { displacement_map[c_v_handle].first + disp[0], displacement_map[c_v_handle].second + disp[1] };
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

					cluster_problem<int>::bound_type objective_boundary = { { -MAX_DISP, -MAX_DISP}, {MAX_DISP, MAX_DISP} };

					pagmo::problem problem{cluster_problem{bound_objective, objective_boundary }};
					pagmo::population pop{ problem, 4 };

					//pagmo::algorithm algo{pagmo::nlopt("neldermead")};
					pagmo::algorithm algo{pagmo::pso{20}};
					pop = algo.evolve(pop);
					std::vector<double> best_value = pop.champion_f();
					std::vector<double> best_arg = pop.champion_x();

					for (auto& c_v_handle : respective_vertices_handles) {						
						displacement_map[c_v_handle].first += best_arg[0];
						displacement_map[c_v_handle].second += best_arg[1];						
						Geometries_with_attributes<Boost_Point_2> c_pt(Boost_Point_2(PointTraits<Point_2>::getX(c_v_handle->point()), PointTraits<Point_2>::getY(c_v_handle->point())));
						c_pt.set_double_attribute("dx", displacement_map[c_v_handle].first);
						c_pt.set_double_attribute("dy", displacement_map[c_v_handle].second);
						pts_gvec.add_geometry(c_pt);
					}
				}
				pts_gvec.to_file(params->temp_dir + "/_" + std::to_string(*distance_val_iter) + "vertices_disp.shp");

			}*/
			
			
			auto aligned_gvec = LinearTopology<EK>::FilteredDisplacedGeovectorFromArrangment(arr, [&](LinearTopology<EK>::Arrangement_2::Vertex_const_iterator viter) {return true; }, agg_displacement_map);
			//aligned_gvec.to_file("C:/DATA_SANDBOX/fault_to_align/data_mo_TIZ110723/data_mo_TIZ110723/example1Bishop/output_dir/aligned_parts.shp");

			GeoVector<Boost_LineString_2> pooled_gvec;
			pool_geovector(aligned_gvec, pooled_gvec, reconstruction_fn);
			pooled_gvec.to_file(params->temp_dir + "/aligned_full.shp");

		}

	}
}