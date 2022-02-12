#include "polygon_decompoe_utils.h"


namespace LxGeo
{
	namespace lxProximityAlign
	{

		SupportPoints decompose_polygons(std::vector<Boost_Polygon_2>& input_polygons, SupportPointsStrategy decompose_strategy) {
			
			SupportPoints out_support_points;

			if (decompose_strategy == SupportPointsStrategy::vertex_and_mid_point) {
				// reserve vectors memory
				out_support_points.polygon_indices.reserve(Constants::MEAN_PTS_PER_POLYGON * 2);
				out_support_points.support_points.reserve(Constants::MEAN_PTS_PER_POLYGON * 2);
				
				size_t c_polygon_idx = 0;
				for (auto& c_polygon : input_polygons) {
					std::list<Boost_Ring_2*> c_polygon_rings;
					c_polygon_rings.push_back(&c_polygon.outer());
					for (Boost_Ring_2& c_inner_ring : c_polygon.inners()) c_polygon_rings.push_back(&c_inner_ring);

					for (auto c_ring : c_polygon_rings) {

						for (size_t c_pt_idx = 0; c_pt_idx < c_ring->size() - 1; ++c_pt_idx) {
							Boost_Point_2& edge_start = c_ring->at(c_pt_idx);
							Boost_Point_2& edge_end = c_ring->at(c_pt_idx+1);
							Boost_Point_2 edge_mid(edge_start.get<0>() + (edge_end.get<0>() - edge_start.get<0>())/2.0f,
								edge_start.get<1>() + (edge_end.get<1>() - edge_start.get<1>())/2.0f );
							out_support_points.polygon_indices.insert(out_support_points.polygon_indices.end(), 3, c_polygon_idx);
							out_support_points.support_points.insert(out_support_points.support_points.end(), { edge_start, edge_mid, edge_end });
						}
					}
					c_polygon_idx++;
				}
				out_support_points.polygon_indices.shrink_to_fit(); out_support_points.support_points.shrink_to_fit();
				return out_support_points;
			}
			else {
				std::cout << "Only vertex_and_mid_point is implemented!" << std::endl;
				throw std::exception("Only vertex_and_mid_point is implemented!");
			}

		}
	}
}