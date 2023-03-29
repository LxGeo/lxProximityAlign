#include "alignment_basic_optim.h"
#include "polygon_decompoe_utils.h"
#include "proximity_triplet.h"
#include "affine_geometry/affine_transformer.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> alignmentIteration(std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster, std::vector<Geometries_with_attributes<Boost_Polygon_2>>& input_polygons) {

			// Support point generation
			SupportPoints c_sup_pts = decompose_polygons(input_polygons, SupportPointsStrategy::constant_walker);

			// spatial to pixel coords
			std::vector<PixelCoords> c_sup_pixels; c_sup_pixels.reserve(c_sup_pts.polygon_indices().size());
			std::transform(c_sup_pts.support_points().begin(), c_sup_pts.support_points().end(),
				std::back_inserter(c_sup_pixels),
				[&ref_raster](auto& pt)->auto{return ref_raster.get_pixel_coords(pt); }
			);

			// proximity triplet reader creation (used to read seperate pixel values from image arrays)
			ProximityTripletLoader PTL(matrices_map["proximity"],
				matrices_map["grad_x"],
				matrices_map["grad_y"]);

			// Read proximity triplet
			std::vector<ProximityTriplet> proximity_triplets; proximity_triplets.reserve(c_sup_pixels.size());
			std::transform(c_sup_pixels.begin(), c_sup_pixels.end(), std::back_inserter(proximity_triplets), [&PTL](auto& px)->ProximityTriplet {return PTL.readTripletAt(px); });

			// convert triplet to disparity (disparity is the proximity multiplied by the sign of the grad at X & Y axis)
			std::vector<SpatialCoords> support_points_disp; support_points_disp.reserve(proximity_triplets.size());
			std::transform(proximity_triplets.begin(), proximity_triplets.end(),
				std::back_inserter(support_points_disp),
				//[&](ProximityTriplet& a)->SpatialCoords { return { a.prox_value * sign(a.grad_y), a.prox_value * sign(a.grad_x) }; }
				ptl_aggregator_function
			);

			////// aggeragte disparites by polygon
			std::vector<SpatialCoords> polygon_disp_values = c_sup_pts.aggregate_points_to_polygon<SpatialCoords, SpatialCoords>(support_points_disp, spatial_coords_median_aggregator);

			//// Align geometries
			//std::vector<Boost_Polygon_2> aligned_polygon = translate_geometries(input_polygons, polygon_disp_values);
			std::vector<Boost_Polygon_2> aligned_polygon; aligned_polygon.reserve(input_polygons.size());
			if (polygon_disp_values.size() == 1) {
				bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(polygon_disp_values[0].xc, polygon_disp_values[0].yc);
				std::transform(input_polygons.begin(), input_polygons.end(), std::back_inserter(aligned_polygon), [&trans_obj](auto& el) {return translate_geometry(el.get_definition(), trans_obj); });
			}
			else {
				auto it1 = input_polygons.begin();
				auto it2 = polygon_disp_values.begin();
				for (; it1 != input_polygons.end() && it2 != polygon_disp_values.end(); ++it1, ++it2)
				{
					bg::strategy::transform::translate_transformer<double, 2, 2> trans_obj(it2->xc, it2->yc);
					aligned_polygon.push_back(
						translate_geometry(it1->get_definition(), trans_obj)
					);
				}
			}

			return aligned_polygon;

		}

	}
}