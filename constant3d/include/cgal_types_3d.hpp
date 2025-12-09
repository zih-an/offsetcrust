#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Circular_kernel_intersections.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/assertions_behaviour.h>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Simple_cartesian.h>


namespace Offset3D {
	//typedef CGAL::Exact_predicates_inexact_constructions_kernel					K;
	typedef CGAL::Exact_predicates_exact_constructions_kernel					K;
	typedef CGAL::Regular_triangulation_vertex_base_3<K>						Vbase;
	typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K, Vbase>		Vb;
	typedef CGAL::Regular_triangulation_cell_base_3<K>							Cb;
#ifdef CGAL_LINKED_WITH_TBB
	typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag>	Tds;
#else
	typedef CGAL::Triangulation_data_structure_3<Vb, Cb>						Tds;
#endif //CGAL_LINKED_WITH_TBB
	typedef CGAL::Regular_triangulation_3<K, Tds>								Regular_triangulation;

	typedef K::Point_3															Point;
	typedef K::Weighted_point_3													Weighted_point;
	typedef K::Vector_3															Vector_3;
	//typedef K::Iso_rectangle_3													Iso_rectangle_3;
	typedef K::Segment_3 														Segment_3;
	typedef K::Ray_3 															Ray_3;
	typedef K::Line_3 															Line_3;
	typedef K::FT  																FT;

	typedef Regular_triangulation::Vertex_handle                                Vertex_handle;
	typedef Regular_triangulation::Cell_handle									Cell_handle;
	typedef Regular_triangulation::Cell_circulator								Cell_circulator;
	typedef Regular_triangulation::Finite_cells_iterator						Finite_cells_iterator;
	typedef Regular_triangulation::Finite_edges_iterator						Edge_iterator;
	typedef Regular_triangulation::Finite_facets_iterator						Finite_facets_iterator;
	typedef Regular_triangulation::All_edges_iterator							All_Edge_iterator;

	//circular geometry
	typedef CGAL::Exact_spherical_kernel_3										Spherical_k;
	typedef Spherical_k::Sphere_3												SSphere_3;
	typedef Spherical_k::Point_3												SPoint_3;

	typedef Spherical_k::Segment_3												SSegment_3;
	typedef Spherical_k::Line_3													SLine_3;
	typedef Spherical_k::Line_arc_3												LineArc3;
	typedef Spherical_k::Circular_arc_3											Sircular_arc_3;
	typedef Spherical_k::Circular_arc_point_3									Sircular_arc_point;


	typedef CGAL::Simple_cartesian<double>										simple_K;
	typedef simple_K::Triangle_3  												Triangle;
	typedef std::list<Triangle>::iterator Iterator;
	typedef CGAL::AABB_triangle_primitive<simple_K, Iterator> Primitive;
	typedef CGAL::AABB_traits<simple_K, Primitive> AABB_triangle_traits;
	typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

	//typedef std::list<Triangle>::iterator 										Iterator;
	//typedef CGAL::AABB_triangle_primitive<K, Iterator> 							Primitive;
	//typedef CGAL::AABB_traits<K, Primitive> 									AABB_triangle_traits;
	//typedef CGAL::AABB_tree<AABB_triangle_traits> 								Tree;

	namespace PMP = CGAL::Polygon_mesh_processing;

}
