#include "alignment_iterative_support_optim.h"
#include "polygon_decompoe_utils.h"
#include "proximity_triplet.h"
#include "affine_geometry/affine_transformer.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> iterative_support_alignment(std::map<std::string, RasterIO>& rasters_map, std::vector<Boost_Polygon_2>& input_polygons) {
			
			size_t N_ITERATION = 1;
			RasterIO& ref_raster = rasters_map["proximity"];
			// proximity triplet reader creation (used to read seperate pixel values from image arrays)
			ProximityTripletLoader PTL(rasters_map["proximity"].raster_data,
				rasters_map["grad_x"].raster_data,
				rasters_map["grad_y"].raster_data);

			// Support point generation
			SupportPoints c_sup_pts = decompose_polygons(input_polygons);

			std::vector<Boost_Point_2> support_points = c_sup_pts.support_points();
			std::vector<SpatialCoords> support_points_disp; support_points_disp.reserve(support_points.size());

			for (size_t c_iteration = 0; c_iteration < N_ITERATION; ++c_iteration) {
				support_points_disp.clear();
				std::cout << "Iteration N: " << c_iteration << std::endl;
				// spatial to pixel coords
				std::vector<PixelCoords> c_sup_pixels; c_sup_pixels.reserve(support_points.size());
				std::transform(support_points.begin(), support_points.end(),
					std::back_inserter(c_sup_pixels),
					[&ref_raster](auto& pt)->auto{return ref_raster.get_pixel_coords(pt); }
				);

				// Read proximity triplet
				std::vector<ProximityTriplet> proximity_triplets; proximity_triplets.reserve(c_sup_pixels.size());
				std::transform(c_sup_pixels.begin(), c_sup_pixels.end(), std::back_inserter(proximity_triplets), [&PTL](auto& px)->ProximityTriplet {return PTL.readTripletAt(px); });

				// convert triplet to disparity (disparity is the proximity multiplied by the sign of the grad at X & Y axis)
				std::transform(proximity_triplets.begin(), proximity_triplets.end(),
					std::back_inserter(support_points_disp),
					ptl_aggregator_function
					//[&](ProximityTriplet& a)->SpatialCoords { return { a.prox_value * sign(a.grad_y) , a.prox_value * sign(a.grad_x) }; }
				);

				//// Align geometries
				std::vector<Boost_Point_2> aligned_support_points = translate_geometries(support_points, support_points_disp);
				support_points = aligned_support_points;

			}
			
			////// aggeragte disparites by polygon
			std::vector<SpatialCoords> polygon_disp_values = c_sup_pts.aggregate_points_to_polygon<SpatialCoords, SpatialCoords>(support_points_disp, spatial_coords_median_aggregator);
			
			//// Align geometries
			std::vector<Boost_Polygon_2> aligned_polygon = translate_geometries(input_polygons, polygon_disp_values);

			return aligned_polygon;

		}

	}
}