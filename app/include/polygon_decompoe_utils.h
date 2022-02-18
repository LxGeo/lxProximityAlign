#pragma once
#include "defs.h"
#include "relationships/composition_struct.h"

namespace LxGeo
{
	namespace lxProximityAlign
	{

		enum SupportPointsStrategy
		{
			vertex_only = 1 << 0,
			vertex_and_mid_point = 1 << 1,
			constant_walker = 1 << 2
		};

		struct SupportPoints : public compositionStrucure<Boost_Point_2> {
			
			//std::vector<size_t> polygon_indices;
			//std::vector<Boost_Point_2> support_points;
			//size_t polygon_count;

			std::vector<size_t>& polygon_indices() { return this->parents_indices; };
			std::vector<Boost_Point_2>& support_points() { return this->children; };
			size_t& polygon_count() { return this->parent_count; };

			template<typename child_proprety, typename parent_proprety>
			std::vector<parent_proprety> aggregate_points_to_polygon(std::vector<child_proprety>& children_properties,
				std::function< parent_proprety(std::list<child_proprety>) > properties_aggregator) { return aggregate_children_to_parent(children_properties, properties_aggregator); }


		};

		SupportPoints decompose_polygons(std::vector<Boost_Polygon_2>& input_polygons,
			SupportPointsStrategy decompose_strategy = SupportPointsStrategy::vertex_and_mid_point);
		
		

	}
}