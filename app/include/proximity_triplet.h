#pragma once
#include "defs.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{
		struct ProximityTriplet {
			double prox_value;
			double grad_x;
			double grad_y;
		};

		class ProximityTripletLoader {

		public:
			ProximityTripletLoader() {};

			ProximityTripletLoader(matrix& _proximity_matrix, matrix& _grad_x_matrix, matrix& _grad_y_matrix) {
				proximity_matrix = _proximity_matrix;
				grad_x_matrix = _grad_x_matrix;
				grad_y_matrix = _grad_y_matrix;
			}

			~ProximityTripletLoader() {};

			ProximityTriplet readTripletAt(PixelCoords& p_c) {
				ProximityTriplet out_triplet = {0,0,0};
				if (p_c.col<0 || p_c.col>=proximity_matrix.cols || p_c.row <0 || p_c.row>=proximity_matrix.rows) {
					std::cout << "coords out of bounds!\n";
					return out_triplet;
				}
				out_triplet.prox_value = proximity_matrix.at<float>(p_c.row, p_c.col);
				out_triplet.grad_x = grad_x_matrix.at<float>(p_c.row, p_c.col);
				out_triplet.grad_y = grad_y_matrix.at<float>(p_c.row, p_c.col);
				return out_triplet;
			}

		public:

			matrix proximity_matrix;
			matrix grad_x_matrix;
			matrix grad_y_matrix;

		};


		auto ptl_aggregator_function = [](ProximityTriplet& ptl)->SpatialCoords {
			//return { ptl.prox_value * sign(ptl.grad_y), ptl.prox_value * sign(ptl.grad_x) };
			if (ptl.grad_x == 0 && ptl.grad_y == 0) return{ 0,0 };
			double grad_x_sq = ptl.grad_x * ptl.grad_x, grad_y_sq = ptl.grad_y * ptl.grad_y;
			SpatialCoords out_disp{
			ptl.prox_value * sign(ptl.grad_y) * grad_y_sq / (grad_x_sq + grad_y_sq),
			ptl.prox_value * sign(ptl.grad_x) * grad_x_sq / (grad_x_sq + grad_y_sq)
			};
			return out_disp;
		};

	}
}