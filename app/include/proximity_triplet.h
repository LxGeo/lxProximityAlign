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
				if (p_c.col<0 || p_c.col>proximity_matrix.cols || p_c.row <0 || p_c.row>proximity_matrix.rows) {
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
	}
}