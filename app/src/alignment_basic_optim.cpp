#include "alignment_basic_optim.h"
#include "polygon_decompoe_utils.h"
#include "proximity_triplet.h"
#include "affine_geometry/affine_transformer.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> alignmentIteration(std::map<std::string, RasterIO>& rasters_map, std::vector<Boost_Polygon_2>& input_polygons) {

			RasterIO& ref_raster = rasters_map["proximity"];
			// Support point generation
			SupportPoints c_sup_pts = decompose_polygons(input_polygons);

			// spatial to pixel coords
			std::vector<PixelCoords> c_sup_pixels; c_sup_pixels.reserve(c_sup_pts.polygon_indices.size());
			std::transform(c_sup_pts.support_points.begin(), c_sup_pts.support_points.end(),
				std::back_inserter(c_sup_pixels),
				[&ref_raster](auto& pt)->auto{return ref_raster.get_pixel_coords(pt); }
			);

			// proximity triplet reader creation (used to read seperate pixel values from image arrays)
			ProximityTripletLoader PTL(rasters_map["proximity"].raster_data,
				rasters_map["grad_x"].raster_data,
				rasters_map["grad_y"].raster_data);

			// Read proximity triplet
			std::vector<ProximityTriplet> proximity_triplets; proximity_triplets.reserve(c_sup_pixels.size());
			std::transform(c_sup_pixels.begin(), c_sup_pixels.end(), std::back_inserter(proximity_triplets), [&PTL](auto& px)->ProximityTriplet {return PTL.readTripletAt(px); });

			// convert triplet to disparity (disparity is the proximity multiplied by the sign of the grad at X & Y axis)
			std::vector<SpatialCoords> support_points_disp; support_points_disp.reserve(proximity_triplets.size());
			std::transform(proximity_triplets.begin(), proximity_triplets.end(),
				std::back_inserter(support_points_disp),
				[&](ProximityTriplet& a)->SpatialCoords { return { a.prox_value * sign(a.grad_x), a.prox_value * sign(a.grad_y) }; }
			);

			////// aggeragte disparites by polygon
			// compute polygon_disparities
			std::map<size_t, std::pair<SpatialCoords, size_t>> polygons_disp_map; // key is polygon_idx & value is pair of (summed support points disparities) & (count of support points)
			for (size_t sup_pt_idx = 0; sup_pt_idx < c_sup_pts.polygon_indices.size(); ++sup_pt_idx) {
				size_t& polygon_idx = c_sup_pts.polygon_indices[sup_pt_idx];
				if (polygons_disp_map.find(polygon_idx) == polygons_disp_map.end()) {
					polygons_disp_map[polygon_idx] = std::make_pair(SpatialCoords({ 0,0 }), 0);
				}
				else {
					polygons_disp_map[polygon_idx].first += support_points_disp[sup_pt_idx];
					polygons_disp_map[polygon_idx].second++;
				}
			}
			// aggregate
			std::vector<SpatialCoords> polygon_disp_values; polygon_disp_values.reserve(polygons_disp_map.size());
			for (size_t poly_idx = 0; poly_idx < input_polygons.size(); ++poly_idx)
				polygon_disp_values.push_back(polygons_disp_map[poly_idx].first / double(polygons_disp_map[poly_idx].second));

			//// Align geometries
			std::vector<Boost_Polygon_2> aligned_polygon = translate_geometries(input_polygons, polygon_disp_values);

			return aligned_polygon;

		}

	}
}