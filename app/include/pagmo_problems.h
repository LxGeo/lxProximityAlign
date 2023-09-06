#pragma once
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/algorithms/pso.hpp>
#include <pagmo/problem.hpp>
#include "defs.h"
#include "parameters.h"
#include "defs_opencv.h"
#include "lightweight/geovector.h"
#include "numcpp/stats.h"
#include "topology/topology_datastructure.h"
#include "graph_weights/spatial_weights.h"
#include "stitching/vector_on_raster_stitcher.h"
#include "design_pattern/defaultmap.h"

namespace LxGeo
{
	using namespace IO_DATA;

	namespace lxProximityAlign
	{

		template <typename Func>
		struct cluster_problem {
			typedef std::pair<std::vector<double>, std::vector<double>> bound_type;
		public:
			cluster_problem() {};
			cluster_problem(Func& _f, bound_type  _bounds) : f(_f), bounds(_bounds) {};
			Func& f;
			bound_type bounds;

			std::vector<double> fitness(const std::vector<double>& dv) const {
				return { f(dv) };
			}
			bound_type get_bounds() const
			{
				return bounds;
			}

		};

		void pagmo_proximity_align_linear(
			std::unordered_map<std::string, matrix>& matrices_map, GeoImage<cv::Mat>& ref_gimg,
			GeoVector<Boost_LineString_2>& input_geovector,
			std::vector<double>& neighbour_distance_band_values,
			std::function<float(numcpp::DetailedStats<float>&)> fitness_from_stats_functor,
			double MAX_DISP,
			std::pair<std::string, std::string> OBJECTIVE_FIELD_NAME_PAIR);

	}
}