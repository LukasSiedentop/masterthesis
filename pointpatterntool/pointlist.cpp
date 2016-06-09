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

pointlist::pointlist(coordinate box, bool periodicity, std::string n) :
		periodic(periodicity), extend(box), name(n) {
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

		list.push_back(new coordinate(0, 0, 0));

		// Tetrahedron
		/*
		 list.push_back(new coordinate(-1, 0, -1));
		 list.push_back(new coordinate(1, 0, -1));
		 list.push_back(new coordinate(0, -1, 1));
		 list.push_back(new coordinate(0, 1, 1));
		 */
		/*
		 list.push_back(new coordinate(0, 2, -1));
		 list.push_back(new coordinate(0, -2, -1));
		 list.push_back(new coordinate(2, 0, 1));
		 list.push_back(new coordinate(-2, 0, 1));
		 */
		/*
		 list.push_back(new coordinate(0, 0, 0));

		 list.push_back(new coordinate(-1, -1, -1));
		 list.push_back(new coordinate(-1, -1, 1));
		 list.push_back(new coordinate(-1, 1, -1));
		 list.push_back(new coordinate(-1, 1, 1));
		 list.push_back(new coordinate(1, -1, -1));
		 list.push_back(new coordinate(1, -1, 1));
		 list.push_back(new coordinate(1, 1, -1));
		 list.push_back(new coordinate(1, 1, 1));
		 */
		extend = coordinate(5, 5, 5);

		break;
	}
	}

	std::cout << name << " point pattern generated." << std::endl;
}

void pointlist::add(double x, double y, double z) {
	coordinate* newPos = new coordinate(x, y, z);

	// ensure that point lies in between boundaries if periodic
	if (periodic) {
		std::vector<coordinate> shifters = extend.getShifters();
		for (std::vector<coordinate>::iterator shifter = shifters.begin();
				shifter != shifters.end(); ++shifter) {
			// as bounding box
			if ((fabs((*shifter)[0] + newPos->x()) < extend.x() / 2)
					&& (fabs((*shifter)[1] + newPos->y()) < extend.y() / 2)
					&& (fabs((*shifter)[2] + newPos->z()) < extend.z() / 2)) {
				(*newPos) += (*shifter);
			}
		}
	}

	// check if point at given position exist
	for (std::vector<coordinate*>::iterator it = list.begin(); it != list.end();
			++it) {
		if (*(*it) == (*newPos)) {
			return;
		}
	}

	// add new point
	list.push_back(newPos);
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

	// example1: http://doc.cgal.org/Manual/3.2/doc_html/cgal_manual/Triangulation_3/Chapter_main.html#Section_22.6
	// example2: http://doc.cgal.org/latest/Geomview/index.html#cite-text-0-0

	nodelist* nlist = new nodelist(periodic, name);

	// initialize geomview for visulaization
	coordinate min = extend / -2, max = extend / 2;
	nlist->setMins(min);
	nlist->setMaxs(max);

	CGAL::Geomview_stream gv(
			CGAL::Bbox_3(min.x(), min.y(), min.z(), max.x(), max.y(), max.z()));
	gv.set_line_width(4);
	// gv.set_trace(true);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	gv.clear();

	if (periodic) {
		// use http://doc.cgal.org/latest/Periodic_3_triangulation_3/index.html
		// calculate delaunay triangulation
		PDT::Iso_cuboid box(min.x(), min.y(), min.z(), max.x(), max.y(),
				max.z());

		PDT PD3d(box);
		std::vector<Point3> points(this->size());

		boost::progress_display show_progress(list.size());

		for (std::vector<coordinate*>::iterator c = list.begin();
				c != list.end(); ++c) {
			points.push_back(Point3((*c)->x(), (*c)->y(), (*c)->z()));
			//std::cout << (*c)->toString() << std::endl;
			//PD3d.insert(Point3((*c)->x(), (*c)->y(), (*c)->z()));
			++show_progress;

			if (list.size() < 100) {
				gv
						<< Sphere3(Point3((*c)->x(), (*c)->y(), (*c)->z()),
								(double) 0.1);
			}

		}
		//PDT PD3d(points.begin(), points.end(), box);
		PD3d.insert(points.begin(), points.end());
		/*
		 gv.set_wired(true);
		 gv << PD3d;
		 gv << CGAL::GREEN;
		 PD3d.draw_dual(gv);
		 gv.set_wired(false);
		 */

		gv << CGAL::RED;
		//PD3d.draw_dual(gv);
		gv << CGAL::GREEN;
		//gv.set_wired(true);
		int i = 0;
		for (PDT::Periodic_tetrahedron_iterator ptit =
				PD3d.periodic_tetrahedra_begin(PDT::UNIQUE_COVER_DOMAIN);
				ptit != PD3d.periodic_tetrahedra_end(PDT::UNIQUE_COVER_DOMAIN);
				++ptit) {

			//PD3d.construct_tetrahedron((*ptit)).vertex(0);
			Point3 centre = CGAL::centroid(
					PD3d.construct_tetrahedron((*ptit)).vertex(0),
					PD3d.construct_tetrahedron((*ptit)).vertex(1),
					PD3d.construct_tetrahedron((*ptit)).vertex(2),
					PD3d.construct_tetrahedron((*ptit)).vertex(3));

			if (i <5) {
				gv << CGAL::GREEN;
				gv << PD3d.construct_tetrahedron((*ptit));
				gv << CGAL::RED;
				gv << Sphere3(centre, 2);
			}

			node* n = nlist->add(centre.x(), centre.y(), centre.z());
			for (int neighIdx = 0; neighIdx < 4; neighIdx++) {
				// exclude infinite cells
				if (!PD3d.is_infinite(ptit.get_cell()->neighbor(neighIdx))) {
					//ptit.get_cell()->neighbor(neighIdx);

					// Neighboiuring cell:
					//ptit.get_cell()->neighbor(neighIdx);

					Point3 neighbour =
							CGAL::centroid(
									ptit.get_cell()->neighbor(neighIdx)->vertex(
											0)->point(),
									ptit.get_cell()->neighbor(neighIdx)->vertex(
											1)->point(),
									ptit.get_cell()->neighbor(neighIdx)->vertex(
											2)->point(),
									ptit.get_cell()->neighbor(neighIdx)->vertex(
											3)->point());
					/*
					 gv << CGAL::BLUE;
					 gv << Sphere3(neighbour, (double) 0.005);
					 */
					node* neigh = nlist->add(neighbour.x(), neighbour.y(),
							neighbour.z());

					// set the neighbourhood
					n->addNeighbour(neigh);
					neigh->addNeighbour(n);
					//}

					if (i ==5) {
						gv << CGAL::BLUE;

						gv << Sphere3(neighbour, 2);
						gv << PD3d.construct_tetrahedron(ptit.get_cell()->neighbor(neighIdx)->vertex(
													0)->point(),
											ptit.get_cell()->neighbor(neighIdx)->vertex(
													1)->point(),
											ptit.get_cell()->neighbor(neighIdx)->vertex(
													2)->point(),
											ptit.get_cell()->neighbor(neighIdx)->vertex(
													3)->point());
					}
				}
			}

			//PD3d.incident_cells();

			i++;
		}
		std::cout << "#cells: " << i << " #cells2: " << PD3d.number_of_cells()
				<< std::endl;
	} else {
		// use http://doc.cgal.org/latest/Triangulation_3/index.html

		// calculate delaunay triangulation
		DT D3d;
		std::vector<Point3> points(this->size());
		for (std::vector<coordinate*>::iterator c = list.begin();
				c != list.end(); ++c) {
			D3d.insert(Point3((*c)->x(), (*c)->y(), (*c)->z()));
			/*
			 if (list.size() < 100) {
			 gv
			 << Sphere3(Point3((*c)->x(), (*c)->y(), (*c)->z()),
			 (double) 0.1);
			 }*/
		}
		/*
		 gv.set_wired(true);
		 gv << D3d;
		 gv << CGAL::GREEN;
		 D3d.draw_dual(gv);
		 gv.set_wired(false);
		 */
		boost::progress_display show_progress(D3d.number_of_cells());

		int i = 0;
		// Iterate over every tetrahedra ("cell" internally)
		for (DT::Finite_cells_iterator cell = D3d.finite_cells_begin();
				cell != D3d.finite_cells_end(); ++cell) {

			// calculate midpoint and add to nodelist
			Point3 centre = CGAL::centroid((*cell).vertex(0)->point(),
					(*cell).vertex(1)->point(), (*cell).vertex(2)->point(),
					(*cell).vertex(3)->point());

			node* n = nlist->add(centre.x(), centre.y(), centre.z());

			/*
			 std::cout << "cell: " << i << std::endl;
			 gv << CGAL::RED;
			 */
			// iterate of all four neighbours, calculate their centroid and set the neighbourhood. TODO: points are being added fourfold...
			for (int neighIdx = 0; neighIdx < 4; neighIdx++) {
				// exclude infinite cells
				if (!D3d.is_infinite((*cell).neighbor(neighIdx))) {

					Point3 neighbour = CGAL::centroid(
							(*cell).neighbor(neighIdx)->vertex(0)->point(),
							(*cell).neighbor(neighIdx)->vertex(1)->point(),
							(*cell).neighbor(neighIdx)->vertex(2)->point(),
							(*cell).neighbor(neighIdx)->vertex(3)->point());
					/*
					 gv << CGAL::BLUE;
					 gv << Sphere3(neighbour, (double) 0.005);
					 */
					node* neigh = nlist->add(neighbour.x(), neighbour.y(),
							neighbour.z());

					// set the neighbourhood
					n->addNeighbour(neigh);
					neigh->addNeighbour(n);
				}
			}

			++show_progress;
			i++;
		}
		//std::cout << "#cells: " << i << " #cells2: " << D3d.number_of_cells() << std::endl;
	}

	t = clock() - t;
	std::cout << "Took " << (((float) t) / CLOCKS_PER_SEC)
			<< "s to generate the neighbours." << std::endl;
	std::cout << "Continue with Enter." << std::endl;
	std::cin.get();
	return nlist;
}
