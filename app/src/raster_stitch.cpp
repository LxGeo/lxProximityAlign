#include "raster_stitch.h"
#include "opencv2/core.hpp"

namespace LxGeo
{
	namespace lxProximityAlign
	{


		double RasterPixelsStitcher::readPolygonsPixels(std::list<Boost_Polygon_2>& resp_polygons, RasterPixelsStitcherStartegy strategy) { 
		

			if (strategy == RasterPixelsStitcherStartegy::contours) {

				double total_obj = 0;
				

				auto line_iter_reader = [](cv::LineIterator& it, matrix& ref_mat)->double {
					double sum = 0;
					const double out_of_extent_ground = 1e5;
					if (*it == nullptr)
						return out_of_extent_ground;
					for (int i = 0; i < it.count; i++, ++it)
						sum += ref_mat.at<float>(it.pos());

					return sum / double(it.count);
				};

				std::list<std::list<cv::Point>> polygons_rings_pixels_coords;
				// rings pixels coords extraction
				for (auto& polygon : resp_polygons) {
					std::list<Boost_Ring_2*> c_polygon_rings;
					c_polygon_rings.push_back(&polygon.outer());
					for (Boost_Ring_2& c_inner_ring : polygon.inners()) c_polygon_rings.push_back(&c_inner_ring);

					for (auto* c_p_ring : c_polygon_rings) {
						std::list<cv::Point> c_ring_pixels;
						for (size_t c_pt_idx = 0; c_pt_idx < c_p_ring->size(); ++c_pt_idx) {
							cv::Point c_pixel_pt;
							ref_raster.get_pixel_coords(c_p_ring->at(c_pt_idx), c_pixel_pt);
							c_ring_pixels.push_back(c_pixel_pt);
						}
						polygons_rings_pixels_coords.push_back(c_ring_pixels);
					}
				}

				for (auto pixel_ring : polygons_rings_pixels_coords) {
					auto pts_iter = pixel_ring.begin();
					for (size_t iter_count = 0; iter_count < pixel_ring.size() - 1; ++iter_count) {
						cv::Point st_pt = (*pts_iter), end_pt = *next(pts_iter);
						if (!ref_raster.pixel_extent.contains(st_pt) || !ref_raster.pixel_extent.contains(end_pt))
							total_obj += 1e5;
						cv::LineIterator it(ref_raster.raster_data, st_pt, end_pt, 8);
						total_obj += line_iter_reader(it, ref_raster.raster_data);
						pts_iter++;
					}
				}

				return total_obj;

			}
			else
				throw std::exception("Only contours strategy is implemented!");


		
		}


	}
}