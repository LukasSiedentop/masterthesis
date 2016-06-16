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

		for (unsigned n = 0; n < 100; n++) {
			// add in a little displacement such that point can not lie directly on the boundary
			list.push_back(
					new coordinate((coordinate(3) * (10)) - (extend) / 2));
		}

		break;
	}
	case 2: { // diamond point pattern ((c) Dirk Ropers)
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

		//list.push_back(new coordinate(0, 0, 0));

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

		// point set large enough so that there are no self edges

		list.push_back(new coordinate(0, 0, 0));

		list.push_back(new coordinate(-1, -1, -1));
		list.push_back(new coordinate(-1, -1, 1));
		list.push_back(new coordinate(-1, 1, -1));
		list.push_back(new coordinate(-1, 1, 1));
		list.push_back(new coordinate(1, -1, -1));
		list.push_back(new coordinate(1, -1, 1));
		list.push_back(new coordinate(1, 1, -1));
		list.push_back(new coordinate(1, 1, 1));

		extend = coordinate(5, 5, 5);

		break;
	}
	}

	std::cout << name << " point pattern generated." << std::endl;
}

coordinate pointlist::centroid(const Tetrahedron tet) const {
	Point3 centr = CGAL::centroid(tet.vertex(0), tet.vertex(1), tet.vertex(2),
			tet.vertex(3));
	return coordinate(centr.x(), centr.y(), centr.z());
}

std::vector<coordinate> pointlist::getNeighbourcentres(Point3 nodeCentre, PDT PD3d) const {

	// **************** debugging ***************
	coordinate min = extend / -2, max = extend / 2;
		CGAL::Geomview_stream gv(
					CGAL::Bbox_3(min.x(), min.y(), min.z(), max.x(), max.y(), max.z()));
			gv.set_line_width(4);
			// gv.set_trace(true);
			gv.set_bg_color(CGAL::Color(200, 200, 0));
			gv.clear();
			// **************** debugging ***************
	std::vector<coordinate> centres;

	 gv << Sphere3(PD3d.locate(nodeCentre)->vertex(0)->point(),0.5);
/*


	// iteration over unique tetrahedra
	PDT::Periodic_tetrahedron pt;
	Tetrahedron tet;
	for (PDT::Periodic_tetrahedron_iterator ptit =
			PD3d.periodic_tetrahedra_begin(PDT::UNIQUE_COVER_DOMAIN);
			ptit != PD3d.periodic_tetrahedra_end(PDT::UNIQUE_COVER_DOMAIN);
			++ptit) {
		pt = (*ptit);
		// Convert the current Periodic_tetrahedron to a Tetrahedron
		tet = PD3d.construct_tetrahedron(pt);
		// construct real tetrahedron

		//tet = PD3d.construct_tetrahedron(pt);

		gv << CGAL::RED;
				gv <<t_bd;
				gv << CGAL::BLUE;
						gv <<tet;

	}
	std::cout << "Continue with Enter." << std::endl;
	std::cin.get();*/
	return centres;
}

Point3 pointlist::getPoint(Periodic_point pt, PDT PD3d) {
	CGAL::Periodic_3_offset_3 offset = pt.second;
	 //PD3d.locate(pt)
	//CGAL::Periodic_3_offset_3 offset =

	for (unsigned int i = 0; i < 3; i++) {

		// shift every offset larger than 2 back
		// 				  n = 0 1 2  3 4 5
		// (((n+1)mod 3)-1) = 0 1 -1 0 1 -1

		offset[i]=(((offset[i]+1)% 3)-1);
	}

	return PD3d.construct_point(pt.first, offset);
}

/*
 std::vector<Point3> pointlist::getPoints(PDT::Cell_iterator cell, PDT PD3d) {
 // inspired by http://cgal-discuss.949826.n4.nabble.com/offset-td950040.html
 // Extract the Periodic_points from the data structure
 Periodic_point pp[4];
 for (int i = 0; i < 4; i++) {
 pp[i] = PD3d.periodic_point(cell, i);
 }

 // Figure out in which directions the cell wraps around
 int off[3];
 for (int d = 0; d < 3; d++) {
 bool does_wrap = false;
 for (int i = 0; i < 4; i++) {
 does_wrap = does_wrap || (pp[i].second[d] != 0);
 }
 off[d] = (does_wrap ? -1 : 0);
 }

 // wrap_offset contains the required translation
 Offset wrap_offset(off[0], off[1], off[2]);

 std::vector<Point3> pts;
 // Apply the translation in wrap_offset to the Periodic_points
 for (int i = 0; i < 4; i++) {
 pp[i] = std::make_pair(pp[i].first, pp[i].second + wrap_offset);
 pts.push_back(PD3d.point(pp[i]));
 }
 return pts;
 }

 bool pointlist::withinOriginalDomain(PDT::Cell_iterator cell, PDT PD3d) {

 for (int i = 0; i < 4; i++) {
 if (!PD3d.periodic_point(cell, i).second.is_null()) {
 return false;
 }
 }
 return true;

 for (int i = 0; i < 4; i++) {
 if (PD3d.periodic_point(cell, i).second.is_null()) {
 return true;
 }
 }
 return false;
 }
 */
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
	//drawGeomview();

	std::cout << "Set the neighbours of the " << name
			<< " point pattern according to the design protocol suggested by Florescu et. al. 2009 (PNAS)."
			<< std::endl;

	// meassure the time
	clock_t t;
	t = clock();

	// example1: http://doc.cgal.org/Manual/3.2/doc_html/cgal_manual/Triangulation_3/Chapter_main.html#Section_22.6
	// example2: http://doc.cgal.org/latest/Geomview/index.html#cite-text-0-0

	coordinate min = extend / -2, max = extend / 2;

	// initialise nodelist
	nodelist* nlist = new nodelist(periodic, name);
	nlist->setMins(min);
	nlist->setMaxs(max);

	// initialize geomview for visulaization
	CGAL::Geomview_stream gv(
			CGAL::Bbox_3(min.x(), min.y(), min.z(), max.x(), max.y(), max.z()));
	gv.set_line_width(4);
	// gv.set_trace(true);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	gv.clear();

	if (periodic) {

		// uses http://doc.cgal.org/latest/Periodic_3_triangulation_3/index.html
		// inspiration: http://cgal-discuss.949826.n4.nabble.com/offset-td950040.html



		PDT::Iso_cuboid domain(min.x(), min.y(), min.z(), max.x(), max.y(),
				max.z());

		PDT PD3d(domain);

		std::cout << "Calculate Delaunay Triangulation..." << std::endl;
		boost::progress_display show_progress(list.size());

		for (std::vector<coordinate*>::iterator c = list.begin();
				c != list.end(); ++c) {
			PD3d.insert(Point3((*c)->x(), (*c)->y(), (*c)->z()));
			// verification
			if (list.size() < 100) {
				gv
						<< Sphere3(Point3((*c)->x(), (*c)->y(), (*c)->z()),
								(double) 0.005);
			}
			++show_progress;
		}

		// iteration over unique tetrahedra that intersects the original domain (centre may lie in between boundaries even if corners of tetrahedron lies outside)
		PDT::Periodic_tetrahedron pt;
		Tetrahedron tet;
		PDT::Cell_handle cell;
		PDT::Cell_handle neighbourCell;

		// estimate duration
		std::cout << "Calculate centroids (nodes) and links to attached tetrahedra (to neighbours)... Takes a lot of time!" << std::endl;
		std::cout << "Iterator estimate: " << PD3d.number_of_cells()*1.5 << std::endl;
		int count = 0;
		show_progress.restart(PD3d.number_of_cells()*1.5);
		for (PDT::Periodic_tetrahedron_iterator ptit =
				PD3d.periodic_tetrahedra_begin(PDT::UNIQUE_COVER_DOMAIN);
				ptit != PD3d.periodic_tetrahedra_end(PDT::UNIQUE_COVER_DOMAIN);
				++ptit) {
			pt = (*ptit);

			// Convert the current Periodic_tetrahedron in T^3 to a Tetrahedron in R^3
			tet = PD3d.construct_tetrahedron(pt);

			// calculate centre
			Point3 nodePoint = CGAL::centroid(tet.vertex(0), tet.vertex(1), tet.vertex(2),
						tet.vertex(3));
			coordinate nodeCentre = coordinate(nodePoint.x(), nodePoint.y(), nodePoint.z());

			// don't evaluate tetrahedron if its centroid is not in the original volume
			if (!nodeCentre.insideAABB(min, max)) {
				continue;
			}

			// add node to nodelist
			node* n = nlist->add(nodeCentre);

			// get cell for the current nodepoint
			cell = PD3d.locate(nodePoint);

			// get neighbouring tetrahedra list: function taking the PD3d and the current tetrahedron. Neighbour is found when three out of four vertices are the same. Calculate centre and add neighbourhoos to nodelist.
			//std::vector<coordinate> neighPoints = getNeighbourcentres(pt, PD3d);
/*
			// **************** debugging *************

			std::cout << "Iterated tetrahedron + centre" << std::endl;
			gv << CGAL::RED;
			gv << Sphere3(tet.vertex(0), 0.01);
			gv << Sphere3(tet.vertex(1), 0.01);
			gv << Sphere3(tet.vertex(2), 0.01);
			gv << Sphere3(tet.vertex(3), 0.01);
			gv << Sphere3(nodePoint, 0.05);

			gv << CGAL::GREEN;
			gv << tet;

			std::cout << "Continue with Enter." << std::endl;
			std::cin.get();

			gv << CGAL::BLUE;

			std::cout << "tetrahedron + centre from cell" << std::endl;
			gv << Sphere3(getPoint(PD3d.periodic_point(PD3d.locate(nodePoint)->vertex(0)), PD3d),0.1);
			gv << Sphere3(getPoint(PD3d.periodic_point(PD3d.locate(nodePoint)->vertex(1)), PD3d),0.1);
			gv << Sphere3(getPoint(PD3d.periodic_point(PD3d.locate(nodePoint)->vertex(2)), PD3d),0.1);
			gv << Sphere3(getPoint(PD3d.periodic_point(PD3d.locate(nodePoint)->vertex(3)), PD3d),0.1);
			//gv << Sphere3(PD3d.point(PD3d.periodic_point(PD3d.locate(pt)->vertex(0))),0.1);

			gv << Sphere3(CGAL::centroid(getPoint(PD3d.periodic_point(PD3d.locate(nodePoint)->vertex(0)), PD3d), getPoint(PD3d.periodic_point(PD3d.locate(nodePoint)->vertex(1)), PD3d), getPoint(PD3d.periodic_point(PD3d.locate(nodePoint)->vertex(2)), PD3d), getPoint(PD3d.periodic_point(PD3d.locate(nodePoint)->vertex(3)), PD3d)),0.5);

			std::cout << "Continue with Enter." << std::endl;
			std::cin.get();
			// **************** debugging *************

*/
			// neighbouriteration
			for (int i = 0; i<4 ; i++) {
				// get current neighbouring cell
				neighbourCell = cell->neighbor(i);

				// construct neighbouring tetrahedron midpoint (by applying the periodicity as in the periodicpointpatterntool)
				Point3 neighbourPoint = CGAL::centroid(getPoint(PD3d.periodic_point(neighbourCell->vertex(0)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(1)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(2)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(3)), PD3d));
				coordinate neighbourCentre = coordinate(neighbourPoint.x(), neighbourPoint.y(), neighbourPoint.z());

				// add the neighbouring node
				node* neighbour = nlist->add(neighbourCentre);

				// TODO: put function addPair(coordinate p1, coordinate p2) in nodelist?
				// set neighbourhood
				n->addNeighbour(neighbour, nlist->getShifter(nodeCentre),
						nlist->getShifter(neighbourCentre));
				neighbour->addNeighbour(n, nlist->getShifter(neighbourCentre),
						nlist->getShifter(nodeCentre));
/*

				// **************** debugging *************
				std::cout << "neighbour " << n << "tetrahedron + centre from cell" << std::endl;
				gv << CGAL::YELLOW;
				gv << Sphere3(getPoint(PD3d.periodic_point(neighbourCell->vertex(0)), PD3d),0.02);
				gv << Sphere3(getPoint(PD3d.periodic_point(neighbourCell->vertex(1)), PD3d),0.02);
				gv << Sphere3(getPoint(PD3d.periodic_point(neighbourCell->vertex(2)), PD3d),0.02);
				gv << Sphere3(getPoint(PD3d.periodic_point(neighbourCell->vertex(3)), PD3d),0.02);

				Tetrahedron neighTet = PD3d.construct_tetrahedron(getPoint(PD3d.periodic_point(neighbourCell->vertex(0)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(1)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(2)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(3)), PD3d));
				gv << neighTet;

				gv << CGAL::GREEN;
				gv << Sphere3(neighbourPoint,0.2);

				std::cout << "Continue with Enter." << std::endl;
				std::cin.get();
				// **************** debugging **************/
			}

			// **************** debugging *************
			//gv.clear();
			// **************** debugging *************
			++show_progress;
			count++;
		}
		std::cout << "Actual size: " << count << std::endl;

		// shift back
		//nlist->shiftList((nlist->getMid()) * (-1));
	} else {
		// use http://doc.cgal.org/latest/Triangulation_3/index.html

		// calculate delaunay triangulation
		DT D3d;
		std::vector<Point3> points(this->size());
		for (std::vector<coordinate*>::iterator c = list.begin();
				c != list.end(); ++c) {
			D3d.insert(Point3((*c)->x(), (*c)->y(), (*c)->z()));

			if (list.size() < 100) {
				gv
						<< Sphere3(Point3((*c)->x(), (*c)->y(), (*c)->z()),
								(double) 0.1);
			}
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

			gv << CGAL::RED;
			gv << Sphere3(centre, (double) 0.005);

			// iterate of all four neighbours, calculate their centroid and set the neighbourhood. TODO: points are being added fourfold...
			for (int neighIdx = 0; neighIdx < 4; neighIdx++) {
				// exclude infinite cells
				if (!D3d.is_infinite((*cell).neighbor(neighIdx))) {

					Point3 neighbour = CGAL::centroid(
							(*cell).neighbor(neighIdx)->vertex(0)->point(),
							(*cell).neighbor(neighIdx)->vertex(1)->point(),
							(*cell).neighbor(neighIdx)->vertex(2)->point(),
							(*cell).neighbor(neighIdx)->vertex(3)->point());

					gv << CGAL::BLUE;
					gv << Sphere3(neighbour, (double) 0.005);

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
	}
	t = clock() - t;
	std::cout << "Took " << (((float) t) / CLOCKS_PER_SEC)
			<< "s to generate the neighbours." << std::endl;
	std::cout << "Continue with Enter." << std::endl;
	std::cin.get();
	return nlist;
}

/*
 void pointlist::drawGeomview() {
 std::vector<std::vector<double> > points;

 for (std::vector<coordinate*>::iterator it = list.begin(); it != list.end();
 ++it) {
 points.push_back(*(*it)->getVector());
 }

 drawPointsGeomview(points, extend / -2, extend / 2);
 }*/
