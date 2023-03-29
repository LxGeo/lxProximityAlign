#include "polygon_decompoe_utils.h"


namespace LxGeo
{
	namespace lxProximityAlign
	{

		SupportPoints decompose_polygons(std::vector<Geometries_with_attributes<Boost_Polygon_2>>& input_polygons, SupportPointsStrategy decompose_strategy) {
			
			SupportPoints out_support_points;

			if (decompose_strategy == SupportPointsStrategy::vertex_and_mid_point) {
				// reserve vectors memory
				out_support_points.polygon_indices().reserve(Constants::MEAN_PTS_PER_POLYGON * 2);
				out_support_points.support_points().reserve(Constants::MEAN_PTS_PER_POLYGON * 2);
				
				size_t c_polygon_idx = 0;
				for (auto& c_polygon_wa : input_polygons) {
					auto& c_polygon = c_polygon_wa.get_definition();
					std::list<Boost_Ring_2*> c_polygon_rings;
					c_polygon_rings.push_back(&c_polygon.outer());
					//for (Boost_Ring_2& c_inner_ring : c_polygon.inners()) c_polygon_rings.push_back(&c_inner_ring);

					for (auto c_ring : c_polygon_rings) {

						for (size_t c_pt_idx = 0; c_pt_idx < c_ring->size() - 1; ++c_pt_idx) {
							Boost_Point_2& edge_start = c_ring->at(c_pt_idx);
							Boost_Point_2& edge_end = c_ring->at(c_pt_idx+1);
							Boost_Point_2 edge_mid(edge_start.get<0>() + (edge_end.get<0>() - edge_start.get<0>())/2.0f,
								edge_start.get<1>() + (edge_end.get<1>() - edge_start.get<1>())/2.0f );
							out_support_points.polygon_indices().insert(out_support_points.polygon_indices().end(), 3, c_polygon_idx);
							out_support_points.support_points().insert(out_support_points.support_points().end(), { edge_start, edge_mid, edge_end });
						}
					}
					c_polygon_idx++;
				}
				out_support_points.polygon_indices().shrink_to_fit(); out_support_points.support_points().shrink_to_fit();
				out_support_points.polygon_count() = c_polygon_idx - 1;
				return out_support_points;
			}
			else if((decompose_strategy == SupportPointsStrategy::constant_walker)){
				double STEP_LENGTH = 5.0;
				// reserve vectors memory
				out_support_points.polygon_indices().reserve(Constants::MEAN_PTS_PER_POLYGON * 10);
				out_support_points.support_points().reserve(Constants::MEAN_PTS_PER_POLYGON * 10);

				size_t c_polygon_idx = 0;
				for (auto& c_polygon_wa : input_polygons) {
					auto& c_polygon = c_polygon_wa.get_definition();
					std::list<Boost_Ring_2*> c_polygon_rings;
					c_polygon_rings.push_back(&c_polygon.outer());
					//for (Boost_Ring_2& c_inner_ring : c_polygon.inners()) c_polygon_rings.push_back(&c_inner_ring);

					for (auto& c_ring : c_polygon_rings) {

						Boost_LineString_2 c_linestring(c_ring->begin(), c_ring->end());
						Boost_MultiPoint_2 densified_ring;
						bg::line_interpolate(c_linestring, STEP_LENGTH, densified_ring);
						out_support_points.polygon_indices().insert(out_support_points.polygon_indices().end(), densified_ring.size(), c_polygon_idx);
						out_support_points.support_points().insert(out_support_points.support_points().end(), densified_ring.begin(), densified_ring.end());
						
					}
					c_polygon_idx++;
				}
				out_support_points.polygon_indices().shrink_to_fit(); out_support_points.support_points().shrink_to_fit();
				out_support_points.polygon_count() = c_polygon_idx - 1;
				return out_support_points;

			}
			else if (decompose_strategy == SupportPointsStrategy::vertex_only) {
				// reserve vectors memory
				out_support_points.polygon_indices().reserve(Constants::MEAN_PTS_PER_POLYGON );
				out_support_points.support_points().reserve(Constants::MEAN_PTS_PER_POLYGON );
				size_t c_polygon_idx = 0;
				for (auto& c_polygon_wa : input_polygons) {
					auto& c_polygon = c_polygon_wa.get_definition();
					std::list<Boost_Ring_2*> c_polygon_rings;
					c_polygon_rings.push_back(&c_polygon.outer());
					//for (Boost_Ring_2& c_inner_ring : c_polygon.inners()) c_polygon_rings.push_back(&c_inner_ring);

					for (auto& c_ring : c_polygon_rings) {
						out_support_points.polygon_indices().insert(out_support_points.polygon_indices().end(), c_ring->size(), c_polygon_idx);
						out_support_points.support_points().insert(out_support_points.support_points().end(), c_ring->begin(), c_ring->end());

					}
					c_polygon_idx++;
				}
				out_support_points.polygon_indices().shrink_to_fit(); out_support_points.support_points().shrink_to_fit();
				out_support_points.polygon_count() = c_polygon_idx - 1;
				return out_support_points;

			}
			else{
				std::cout << "Only vertex_only & vertex_and_mid_point & constant_walker are implemented!" << std::endl;
				throw std::exception("Only vertex_only & vertex_and_mid_point & constant_walker are implemented!");
			}

		}
	}
}