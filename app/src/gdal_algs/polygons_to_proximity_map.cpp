#include "gdal_algs/polygons_to_proximity_map.h"

namespace LxGeo
{
    namespace lxProximityAlign
    {

		using namespace LxGeo::GeometryFactoryShared;

        void polygons2proximity(IO_DATA::PolygonsShapfileIO& input_shapefile, std::string& out_proximity_map_path, OGREnvelope* out_extents,
			double raster_px_size, double raster_py_size, ProximityMapStrategy proximity_strategy) {

			std::string rasterized_out_path = (boost::filesystem::path(params->temp_dir)/"rasterized.tif").string();
			
			double constant_walker_step = 10.0;

            // transform geometries depending on strategy
            std::list<OGRGeometryH> geometries_to_burn;
			input_shapefile.vector_layer->ResetReading();
			tqdm bar;
			for (size_t j = 0; j < input_shapefile.feature_count; ++j) {
				bar.progress(j, input_shapefile.feature_count);
				OGRFeature* feat = input_shapefile.vector_layer->GetNextFeature();
				if (feat == NULL) continue;

				OGRGeometry* geom = feat->GetGeometryRef();

				// Assumes the shapefile only contains OGRPolygons
				if (OGRPolygon* P = dynamic_cast<OGRPolygon*>(geom)) {
					switch (proximity_strategy) {
						case ProximityMapStrategy::vertex_only:
						{
							std::list<OGRLinearRing*> rings;
							rings.push_back(P->getExteriorRing()); for (size_t in_r_idx = 0; in_r_idx < P->getNumInteriorRings(); ++in_r_idx) rings.push_back(P->getInteriorRing(in_r_idx));
							for (auto c_ring : rings) {
								OGRLineString* ring_pts = c_ring->toLineString();
								for (size_t c_r_pt_idx = 0; c_r_pt_idx < c_ring->getNumPoints(); ++c_r_pt_idx) {
									OGRPoint c_pt;
									c_ring->getPoint(c_r_pt_idx, &c_pt);
									geometries_to_burn.push_back(c_pt.clone());
								}
							}
							break;
						}
						case ProximityMapStrategy::contours:
						{
							std::list<OGRLinearRing> rings;
							geometries_to_burn.push_back(P->getExteriorRing()->toLineString());
							for (size_t in_r_idx = 0; in_r_idx < P->getNumInteriorRings(); ++in_r_idx) geometries_to_burn.push_back(P->getInteriorRing(in_r_idx)->toLineString());
							break;
						}
						case ProximityMapStrategy::filled_polygon:
						{
							geometries_to_burn.push_back(P);
							break;
						}
						case ProximityMapStrategy::constant_walker:
						{
							std::list<OGRLinearRing*> rings;
							rings.push_back(P->getExteriorRing()); for (size_t in_r_idx = 0; in_r_idx < P->getNumInteriorRings(); ++in_r_idx) rings.push_back(P->getInteriorRing(in_r_idx));
							for (auto c_ring : rings) {
								OGRLineString* ring_pts = c_ring->toLineString();
								ring_pts->segmentize(constant_walker_step);
								for (size_t c_r_pt_idx = 0; c_r_pt_idx < c_ring->getNumPoints(); ++c_r_pt_idx) {
									OGRPoint c_pt;
									c_ring->getPoint(c_r_pt_idx, &c_pt);
									geometries_to_burn.push_back(c_pt.clone());
								}
							}
							break;
						}
						case ProximityMapStrategy::skeleton:
						{
							geometries_to_burn.push_back(BuildMultiLine(P)->clone());
							break;
						}
					}
				}
			}     
			bar.finish();

			std::string aux_shapefile_path = (boost::filesystem::path(params->temp_dir) / "aux_shapefile.shp").string();
			switch (proximity_strategy) {
				case ProximityMapStrategy::vertex_only:
				case ProximityMapStrategy::constant_walker:
				{					
					IO_DATA::PointsShapfileIO aux_shapefile(aux_shapefile_path, input_shapefile.spatial_refrence);
					aux_shapefile.write_geometries(geometries_to_burn);
					break;
				}
				case ProximityMapStrategy::contours:
				{
					IO_DATA::LineStringShapfileIO aux_shapefile(aux_shapefile_path, input_shapefile.spatial_refrence);
					aux_shapefile.write_geometries(geometries_to_burn);
					break;
				}
				case ProximityMapStrategy::filled_polygon:
				{
					IO_DATA::PolygonsShapfileIO aux_shapefile(aux_shapefile_path, input_shapefile.spatial_refrence);
					aux_shapefile.write_geometries(geometries_to_burn);
					break;
				}
				case ProximityMapStrategy::skeleton:
				{
					IO_DATA::LineStringShapfileIO aux_shapefile(aux_shapefile_path, input_shapefile.spatial_refrence);
					aux_shapefile.write_geometries(geometries_to_burn);
					break;
				}
			}

			rasterize_shapefile(rasterized_out_path, aux_shapefile_path, out_extents, raster_px_size, raster_py_size);

			std::string proximity_out_path = (boost::filesystem::path(params->temp_dir) / "proximity.tif").string();
			transform_to_proximity(rasterized_out_path, proximity_out_path);


        }

		void linestrings2proximity(IO_DATA::LineStringShapfileIO& input_shapefile, std::string& out_proximity_map_path, OGREnvelope* out_extents,
			double raster_px_size, double raster_py_size, ProximityMapStrategy proximity_strategy) {

			std::string rasterized_out_path = (boost::filesystem::path(params->temp_dir) / "rasterized.tif").string();
			const char* shapefilePath = input_shapefile.vector_dataset->GetDescription();
			rasterize_shapefile(rasterized_out_path, std::string(shapefilePath), out_extents, raster_px_size, raster_py_size);

			std::string proximity_out_path = (boost::filesystem::path(params->temp_dir) / "proximity.tif").string();
			transform_to_proximity(rasterized_out_path, proximity_out_path);
		}
    }
}