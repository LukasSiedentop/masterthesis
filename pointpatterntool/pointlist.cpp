/*
 * pointlist.cpp
 *
 *  Created on: 07.06.2016
 *      Author: Lukas Siedentop
 *
 *      Class that describes a simple collection of points from which a decorated pattern (nodelist) can be calculated.
 */

#include "pointlist.hpp"

pointlist::pointlist() :
		periodic(0), extend(coordinate()) {
}

pointlist::pointlist(coordinate Box, bool periodicity, std::string n) :
		periodic(periodicity), extend(Box), name(n) {
}

pointlist::pointlist(int pattern, bool periodicity) :
		periodic(periodicity), extend(coordinate(10, 10, 10)) {

	switch (pattern) {
	case 1: { // poisson pattern
		std::cout
				<< "Generate Poisson pattern with 1000 points in 10^3 cubicle..."
				<< std::endl;

		name = "Poisson";

		// seed for random numbers generator
		srand(time(NULL));

		for (unsigned n = 0; n < 1000; n++) {
			list.push_back(new coordinate((coordinate(3) * 10) - extend / 2));
		}

		break;
	}
	case 2: { // diamond point pattern ((c) Dirk)
		std::cout
				<< "Generate Diamond point pattern with 1000 points in 10^3 cubicle..."
				<< std::endl;

		name = "Diamond";

		// TODO: generate unit cell and then extend to desired size

		int x = -5, y = -5, z = -5;
		int numx = 5, numy = 5, numz = 5;

		while (numz - z > 0) {
			y = -5;
			while (numy - y > 0) {
				x = -5;
				while (numx - x > 0) {
					if (not ((x + y + z) % 2)) {
						// first point: (0.25,0.25,0.25)
						list.push_back(
								new coordinate(
										(coordinate(0.25, 0.25, 0.25)
												+ coordinate(x, y, z))));
						// second point: (0.75,0.75,0.75)
						list.push_back(
								new coordinate(
										(coordinate(0.75, 0.75, 0.75)
												+ coordinate(x, y, z))));
					}
					x++;
				}
				y++;
			}
			z++;
		}

		break;
	}

	case 3: { // test pattern (tetrahedron,...)
		std::cout
				<< "Generate Test point pattern with 9 points in 10^3 cubicle..."
				<< std::endl;

		name = "Test";
		/*
		 // Tetrahedron
		 list.push_back(new node(this, (coordinate(-1, 0, -1))));
		 list.push_back(new node(this, (coordinate(1, 0, -1))));
		 list.push_back(new node(this, (coordinate(0, -1, 1))));
		 list.push_back(new node(this, (coordinate(0, 1, 1))));
		 */

		list.push_back(new coordinate(0, 0, 0));

		list.push_back(new coordinate(-1, -1, -1));
		list.push_back(new coordinate(-1, -1, 1));
		list.push_back(new coordinate(-1, 1, -1));
		list.push_back(new coordinate(-1, 1, 1));
		list.push_back(new coordinate(1, -1, -1));
		list.push_back(new coordinate(1, -1, 1));
		list.push_back(new coordinate(1, 1, -1));
		list.push_back(new coordinate(1, 1, 1));

		extend = coordinate(3, 3, 3);

		break;
	}
	}

	std::cout << name << " point pattern generated." << std::endl;
}

void pointlist::add(double x, double y, double z) {
	// check if point at given position exist
	for (std::vector<coordinate*>::iterator it = list.begin(); it != list.end();
			++it) {
		if (periodic) {
			coordinate pos = *(*it);
			if ((std::min(fabs(x - pos[0]),
					(fabs(fabs(x - pos[0]) - extend.x()))) < 0.000001)
					&& (std::min(fabs(y - pos[1]),
							(fabs(fabs(y - pos[1]) - extend.y()))) < 0.000001)
					&& (std::min(fabs(z - pos[2]),
							(fabs(fabs(z - pos[2]) - extend.z()))) < 0.000001)) {
				return;
			}
		}
		if (!periodic && (*(*it) == coordinate(x, y, z))) {
			return;
		}
	}

	coordinate newPos = coordinate(x, y, z);

	// ensure that point lies in between boundaries if periodic
	if (periodic) {
		std::vector<coordinate> shifters = extend.getShifters();
		for (std::vector<coordinate>::iterator shifter = shifters.begin();
				shifter != shifters.end(); ++shifter) {
			// as bounding box
			if ((fabs((*shifter)[0] + newPos[0]) < extend.x() / 2)
					&& (fabs((*shifter)[1] + newPos[1]) < extend.y() / 2)
					&& (fabs((*shifter)[2] + newPos[2]) < extend.z() / 2)) {
				newPos += (*shifter);
			}
		}
	}

	// build and add new point
	coordinate* c = new coordinate(x, y, z);
	list.push_back(c);
}

void pointlist::setBox(coordinate box) {
	extend = box;
}

unsigned int pointlist::size() const {
	return list.size();
}

void pointlist::deleteEntries() {
	for (std::vector<coordinate*>::iterator c = list.begin(); c != list.end();
			++c) {
		delete (*c);
	}
}

void pointlist::scaleList(double a) {
	// correct extend box
	extend *= a;

	// scale entries
	for (std::vector<coordinate*>::iterator it = list.begin(); it != list.end();
			++it) {
		*(*it) *= a;
	}
}

nodelist* pointlist::decorate() {
	std::cout << "Set the neighbours of the " << name
			<< " point pattern according to the design protocol suggested by Florescu et. al. 2009 (PNAS)."
			<< std::endl;

	// meassure the time
	clock_t t;
	t = clock();
	// TODO:
	// Delaunay triangulation via CGAL
	// connect neighbouring tetrahedron centre of gravities.

	// example1: http://doc.cgal.org/Manual/3.2/doc_html/cgal_manual/Triangulation_3/Chapter_main.html#Section_22.6
	// example2: http://doc.cgal.org/latest/Geomview/index.html#cite-text-0-0

	nodelist* nlist = new nodelist(periodic, name);

	// initialize geomview for visulaization
	coordinate min = extend / -2, max = extend / 2;

	std::cout << "min: " <<  min.toString() << " max: " << max.toString() << std::endl;

	CGAL::Geomview_stream gv(
			CGAL::Bbox_3(min.x(), min.y(), min.z(), max.x(), max.y(), max.z()));
	gv.clear();

	for (std::vector<coordinate*>::iterator c = list.begin(); c != list.end();
			++c) {
		std::vector<double> pt = *(*c)->getVector();

		if (this->size() < 100) {
			// draw point pattern
			gv << Sphere3(Point3((*c)->x(), (*c)->y(), (*c)->z()), (double) 0.01);
		}

	}

	if (periodic) {
		// use http://doc.cgal.org/latest/Periodic_3_triangulation_3/index.html
		// calculate delaunay triangulation
		PDT::Iso_cuboid box(min.x(), min.y(), min.z(), max.x(), max.y(),
				max.z());
		std::vector<Point3> points(this->size());
		PDT PD3d(box);
		std::cout << "Start to add points to the triangulation..." << std::endl;
		int i =1;
		for (std::vector<coordinate*>::iterator c = list.begin();
				c != list.end(); ++c) {
			PD3d.insert(Point3((*c)->x(), (*c)->y(), (*c)->z()));
			i++;
			//points.push_back();
		}

		std::cout << "Periodic triangulation built!" << std::endl;

		// Iterate over every tetrahedron ("cell" internally)
		i =0;
		for (CGAL::Periodic_3_triangulation_3<GT>::Cell_iterator cell = PD3d.cells_begin();
				cell != PD3d.cells_end(); ++cell) {
			i++;
			// iterates also periodic continued pattern...
			// TODO: calculate midpoint, add to nodelist with according neighbours
			//(*cell).
		}
		std::cout << "Tetrahedrons: " << i << " or " << PD3d.number_of_cells() << std::endl;

	} else {
		// use http://doc.cgal.org/latest/Triangulation_3/index.html

		// calculate delaunay triangulation
		DT D3d;
		std::vector<Point3> points(this->size());
		for (std::vector<coordinate*>::iterator c = list.begin();
				c != list.end(); ++c) {
			D3d.insert(Point3((*c)->x(), (*c)->y(), (*c)->z()));
		}

		std::cout << "Non-periodic triangulation built!" << std::endl;

		// Iterate over every tetrahedron ("cell" internally)
		int i =0;
		for (CGAL::Triangulation_3<K>::Cell_iterator cell = D3d.cells_begin();
				cell != D3d.cells_end(); ++cell) {
			i++;
			// TODO: calculate midpoint, add to nodelist with according neighbours
			//(*cell).
		}
		std::cout << "Tetrahedrons: " << i << " or " << D3d.number_of_cells() << std::endl;
		/*
		 std::cout << "centroid point pattern: "
		 << CGAL::centroid(points.begin(), points.end()) << std::endl;
		 */
		std::cout << "Drawing 3D Delaunay triangulation in wired mode."
				<< std::endl;
		//gv << CGAL::BLUE;
		gv.set_wired(true);
		gv << D3d;

		// voronoi pattern
		//gv << CGAL::RED;
		//gv.set_line_width(4);
		//D3d.draw_dual(gv);

	}

	std::cout << "Weiter mit Enter." << std::endl;
	std::cin.get();
	//boost::progress_display show_progress(1);
	//++show_progress;

	t = clock() - t;
	std::cout << "Took " << (((float) t) / CLOCKS_PER_SEC)
			<< "s to generate the neighbours." << std::endl;

	return nlist;
}
