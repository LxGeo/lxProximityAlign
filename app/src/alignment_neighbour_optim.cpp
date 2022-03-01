#include "defs.h"
#include "parameters.h"
#include "io_shapefile.h"
#include "alignment_neighbour_optim.h"
#include "polygon_decompoe_utils.h"
#include "proximity_triplet.h"
#include "affine_geometry/affine_transformer.h"
#include "graph_weights/polygons_spatial_weights.h"
#include "graph_weights/spatial_weights.h"
#include "relationships/composition_struct.h"

namespace LxGeo
{

	using namespace IO_DATA;
	namespace lxProximityAlign
	{

		std::vector<Boost_Polygon_2> alignmentNeighbour(std::map<std::string, matrix>& matrices_map, RasterIO& ref_raster, std::vector<Boost_Polygon_2>& input_polygons) {

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
				ptl_aggregator_function
			);

			////// aggeragte disparites by polygon
			std::vector<SpatialCoords> polygon_disp_values = c_sup_pts.aggregate_points_to_polygon<SpatialCoords, SpatialCoords>(support_points_disp, spatial_coords_mean_aggregator);
			
			// maybe add linetracer to weight polygons

			// Spatial weights
			PolygonSpatialWeights PSW = PolygonSpatialWeights(input_polygons);
			WeightsDistanceBandParams wdbp = { 5, false, -1, [](double x)->double { return 1.0 / (1.0 + x); } };
			PSW.fill_distance_band_graph(wdbp);
			PSW.run_labeling();

			std::vector<LineString_with_attributes> edges_line_strings = PSW.export_edge_graph_as_LSwithAttr();
			std::string component_shapefile_path = (boost::filesystem::path(params->temp_dir) / "components.shp").string();
			LineStringShapfileIO edges_out_shapefile = LineStringShapfileIO(component_shapefile_path, ref_raster.spatial_refrence);
			edges_out_shapefile.write_linestring_shapefile(edges_line_strings);

			// aggregate diparites by block
			compositionStrucure<Boost_Polygon_2> block_to_poly_relation; 
			block_to_poly_relation.parents_indices = PSW.component_labels;
			block_to_poly_relation.parent_count = PSW.n_components;
			std::vector<SpatialCoords> block_disp_values= block_to_poly_relation.aggregate_children_to_parent<SpatialCoords, SpatialCoords>(polygon_disp_values, spatial_coords_mean_aggregator);

			for (size_t poly_idx = 0; poly_idx < polygon_disp_values.size(); ++poly_idx) {
				size_t respective_block_idx = PSW.component_labels[poly_idx];
				polygon_disp_values[poly_idx] = block_disp_values[respective_block_idx];
			}

			//// Align geometries
			std::vector<Boost_Polygon_2> aligned_polygon = translate_geometries(input_polygons, polygon_disp_values);

			return aligned_polygon;
		}

	}
}