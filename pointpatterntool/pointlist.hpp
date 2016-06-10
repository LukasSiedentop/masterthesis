/*
 * pointlist.hpp
 *
 *  Created on: 07.06.2016
 *      Author: Lukas Siedentop
 *
 *      Class that describes a simple collection of points from which a decorated pattern (nodelist) can be calculated.
 *      Also operations and statistics based on the pure pattern should be calculated here.
 */

#ifndef POINTLIST_HPP_
#define POINTLIST_HPP_

#include "coordinate.hpp"
#include "nodelist.hpp"
//#include "functions.hpp"

// used for design protocol
#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/Delaunay_triangulation_3.h> // calculate Delaunay triangulation from set  of points

#include <CGAL/IO/Geomview_stream.h> // geomview
#include <CGAL/IO/Triangulation_geomview_ostream_3.h> // geomview ability to draw 3d-triangulations

#include <CGAL/centroid.h> // calculate centroid from set of points

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>

#include <CGAL/Periodic_3_triangulation_ds_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Periodic_3_triangulation_ds_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>

// CGAL
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Cartesian<double>  K;

// 3D triangulation data structure
typedef CGAL::Triangulation_data_structure_3<>	TDS;

typedef TDS::Cell_handle	Cell_handle;

// non-periodic delaunay triangulation
typedef CGAL::Delaunay_triangulation_3<K>   DT;
typedef K::Point_3 Point3;
typedef K::Sphere_3 Sphere3;
//typedef K::Segment_3 Segment3;

// periodic delaunay triangulation
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K> PT;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<PT> PDT;
//typedef GT::Sphere_3 Pshere3;

typedef PDT::Periodic_point         Periodic_point;
typedef PDT::Offset                 Offset;

class pointlist {
private:
	bool periodic;

	// sidelengths
	coordinate extend;

	std::string name;

	std::vector<class coordinate*> list;

	// returns the points of a given periodic data structure cell, periodically shifted such that the original cell is centred
	std::vector<Point3> getPoints(PDT::Cell_iterator cell, PDT PD3d);
public:
	pointlist();
	pointlist(coordinate box, bool periodicity, std::string n);
	// constructs a pattern (density of points 1, within 10^3 cubicle) with: pattern=1 - random points, pattern=2 - points arranged in a diamond lattice, pattern=3 - explicit testpoints.
	pointlist(int pattern, bool periodicity);

	void setBox(coordinate box);

	unsigned int size() const;

	// Adds a point at the given coordinate if none exists there
	void add(double x, double y, double z);
	void deleteEntries();

	void scaleList(double a);

	// returns the decorated pattern, constructed with the design protocol by Florescu et. al. 2009 (PNAS)
	nodelist* decorate();

	/* TODO:
	 * Structure factor
	 * Hyperuniformity
	 * plot geomview & gnuplot
	 */
};

#endif /* POINTLIST_HPP_ */
