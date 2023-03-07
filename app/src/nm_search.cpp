#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "nm_search.h"
#include "optim.hpp"
#include "numcpp/line_space.h"

#include "graph_weights/polygons_spatial_weights.h"
#include "graph_weights/spatial_weights.h"
#include "affine_geometry/affine_transformer.h"
#include "raster_stitch.h"
#include "io_shapefile.h"
#include "parameters.h"


namespace LxGeo
{
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> nm_proximity_align(std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster, std::vector<Boost_Polygon_2>& input_polygons) {

			std::vector<Boost_Polygon_2> out_geometries(input_polygons.begin(), input_polygons.end());
			std::vector<double> geometry_traversal_weights(input_polygons.size(), DBL_MAX);

			std::vector<double> neighbour_distance_band_values = {0.01, 1, 5, 10, 20};//numcpp::linspace(0.0, 10, 9);
			PolygonSpatialWeights PSW = PolygonSpatialWeights(input_polygons);
			WeightsDistanceBandParams wdbp = { 20, false, -1, [](double x)->double { return x; } };
			PSW.fill_distance_band_graph(wdbp);

			RasterPixelsStitcher RPR(ref_raster);

			auto objective_fn = [&RPR](const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)->double {
				std::list<SpatialCoords> c_displacement = { { vals_inp(0), vals_inp(1) } };
				auto resp_polys = static_cast<std::list<Boost_Polygon_2>*>(opt_data);
				std::list<Boost_Polygon_2> translated_resp_polys= translate_geometries(*resp_polys, c_displacement);
				double obj_val = RPR.readPolygonsPixels(translated_resp_polys, RasterPixelsStitcherStartegy::contours);
				return obj_val;
			};

			
			double MAX_DISP=params->max_disparity;
			for (auto& distance_val_iter = neighbour_distance_band_values.rbegin(); distance_val_iter != neighbour_distance_band_values.rend(); ++distance_val_iter) {

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
					std::list<Boost_Polygon_2> respective_polygons; for (auto poly_idx : components_polygons_map[comp_idx]) respective_polygons.push_back(out_geometries[poly_idx]);
					
					optim::algo_settings_t c_settings;
					c_settings.iter_max = 50;
					c_settings.vals_bound = true; c_settings.lower_bounds = -MAX_DISP * Eigen::VectorXd::Ones(2); c_settings.upper_bounds = MAX_DISP * Eigen::VectorXd::Ones(2);
					//c_settings.print_level = 3;
					optim::Mat_t simplex_points(3, 2);
					simplex_points.row(0) << MAX_DISP * ((double)rand() / RAND_MAX - 0.5), MAX_DISP* ((double)rand() / RAND_MAX - 0.5); // 1.5, 1.5;
					simplex_points.row(1) << MAX_DISP * ((double)rand() / RAND_MAX - 0.5), MAX_DISP* ((double)rand() / RAND_MAX - 0.5); //0.1, 0.15;
					simplex_points.row(2) << MAX_DISP * ((double)rand() / RAND_MAX - 0.5), MAX_DISP* ((double)rand() / RAND_MAX - 0.5); //-1.0, 1.0;
					c_settings.nm_settings.custom_initial_simplex = true; c_settings.nm_settings.initial_simplex_points = simplex_points;

					bool success = optim::nm(init_vals, objective_fn, &respective_polygons, c_settings);
					if (!success) { std::cout << "optimization failed!" << std::endl; continue; }

					for (auto poly_idx : components_polygons_map[comp_idx]) {
						bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(init_vals(0), init_vals(1));
						out_geometries[poly_idx] = translate_geometry(out_geometries[poly_idx], trans_obj);
					}
								
				}
				bar.finish();
				//break;
				/*
				PolygonsShapfileIO aligned_out_shapefile = PolygonsShapfileIO(params->temp_dir + "/dist_"+std::to_string(*distance_val_iter)+".shp", nullptr);
				auto polygons_with_attrs = transform_to_geom_with_attr<Boost_Polygon_2>(out_geometries);
				std::cout << "Writing outfile!" << std::endl;
				aligned_out_shapefile.write_shapefile(polygons_with_attrs);*/

			}
			return out_geometries;

		}

		std::vector<Boost_Polygon_2> nm_proximity_align_1d(
			std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster,
			std::vector<Boost_Polygon_2>& input_polygons,
			std::pair<double, double>& r2r_constants
		) {

			std::vector<Boost_Polygon_2> out_geometries(input_polygons.begin(), input_polygons.end());

			std::vector<double> neighbour_distance_band_values = { 0.01, 1, 5, 10, 20 };
			PolygonSpatialWeights PSW = PolygonSpatialWeights(input_polygons);
			WeightsDistanceBandParams wdbp = { 20, false, -1, [](double x)->double { return x; } };
			PSW.fill_distance_band_graph(wdbp);

			RasterPixelsStitcher RPR(ref_raster);

			auto objective_fn = [&RPR, &r2r_constants](const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)->double {
				auto height = vals_inp(0);
				std::list<SpatialCoords> c_displacement = { { r2r_constants.first* height, r2r_constants.second * height } };
				auto resp_polys = static_cast<std::list<Boost_Polygon_2>*>(opt_data);
				std::list<Boost_Polygon_2> translated_resp_polys = translate_geometries(*resp_polys, c_displacement);
				double obj_val = RPR.readPolygonsPixels(translated_resp_polys, RasterPixelsStitcherStartegy::contours);
				return obj_val;
			};

			double MAX_DISP = params->max_disparity;
			for (auto& distance_val_iter = neighbour_distance_band_values.rbegin(); distance_val_iter != neighbour_distance_band_values.rend(); ++distance_val_iter) {

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
				for (size_t comp_idx = 0; comp_idx < PSW.n_components; ++comp_idx) {
					bar.progress(comp_idx, PSW.n_components);
					auto& init_vals = components_transforms_map[comp_idx];
					std::list<Boost_Polygon_2> respective_polygons; for (auto poly_idx : components_polygons_map[comp_idx]) respective_polygons.push_back(out_geometries[poly_idx]);

					optim::algo_settings_t c_settings;
					c_settings.iter_max = 50;
					c_settings.vals_bound = true; c_settings.lower_bounds = Eigen::VectorXd::Zero(1); c_settings.upper_bounds = MAX_DISP * Eigen::VectorXd::Ones(1);
					//c_settings.print_level = 3;
					optim::Mat_t simplex_points(2, 1);
					simplex_points.row(0) << MAX_DISP * ((double)rand() / RAND_MAX ); // 1.5, 1.5;
					simplex_points.row(1) << MAX_DISP * ((double)rand() / RAND_MAX ); //0.1, 0.15;
					//simplex_points.row(2) << MAX_DISP * ((double)rand() / RAND_MAX - 0.5); //-1.0, 1.0;
					c_settings.nm_settings.custom_initial_simplex = true; c_settings.nm_settings.initial_simplex_points = simplex_points;

					bool success = optim::nm(init_vals, objective_fn, &respective_polygons, c_settings);
					if (!success) { std::cout << "optimization failed!" << std::endl; continue; }

					for (auto poly_idx : components_polygons_map[comp_idx]) {
						bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(r2r_constants.first * init_vals(0), r2r_constants.second*init_vals(0));
						out_geometries[poly_idx] = translate_geometry(out_geometries[poly_idx], trans_obj);
					}

				}
				bar.finish();
			}
			return out_geometries;

		}

	}
}