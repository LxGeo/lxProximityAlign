#pragma once
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "defs.h"
#include "parameters.h"
#include "io_raster.h"
#include "geometries_with_attributes/geometries_with_attributes.h"
#include "lightweight/geovector.h"
#include "numcpp/stats.h"
#include "optim.hpp"
#include "numcpp/line_space.h"
#include "graph_weights/spatial_weights.h"
#include "affine_geometry/affine_transformer.h"
#include "stitching/vector_on_raster_stitcher.h"
#include "io_shapefile.h"
#include "optim_algo/min_search_golden_section.h"
#include "optim_algo/particle_swarm.h"
#include "topology/topology_datastructure.h"
#include "design_pattern/cache.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		void nm_proximity_align(std::unordered_map<std::string, matrix>& matrices_map, GeoImage<cv::Mat>& ref_gimg,
			GeoVector<Boost_Polygon_2>& input_geovector,
			std::vector<double>& neighbour_distance_band_values,
			std::function<float(numcpp::DetailedStats<float>&)> fitness_from_stats_functor,
			double MAX_DISP,
			std::pair<std::string, std::string> OBJECTIVE_FIELD_NAME_PAIR = {"DISPARIT_X", "DISPARIT_Y"});

		void nm_proximity_align_1d(
			std::unordered_map<std::string, matrix>& matrices_map, GeoImage<cv::Mat>& ref_raster,
			GeoVector<Boost_Polygon_2>& input_polygons,
			std::vector<double>& neighbour_distance_band_values,
			std::pair<double, double>& r2r_constants,
			std::function<float(numcpp::DetailedStats<float>&)> fitness_from_stats_functor,
			double MAX_DISP = 100,
			std::string DISP_COLUMN_NAME = "DISP"
		);

		void nm_proximity_align_linear(
			std::unordered_map<std::string, matrix>& matrices_map, GeoImage<cv::Mat>& ref_gimg,
			GeoVector<Boost_LineString_2>& input_geovector,
			std::vector<double>& neighbour_distance_band_values,
			std::function<float(numcpp::DetailedStats<float>&)> fitness_from_stats_functor,
			double MAX_DISP,
			std::pair<std::string, std::string> OBJECTIVE_FIELD_NAME_PAIR);

	}
}