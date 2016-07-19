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
		periodic(0), min(coordinate()), max(coordinate()) {
}

pointlist::pointlist(coordinate boxMin, coordinate boxMax, bool periodicity,
		std::string n) :
		periodic(periodicity), min(boxMin), max(boxMax), name(n) {
}

pointlist::pointlist(int pattern, bool periodicity) :
		periodic(periodicity), min(coordinate(-5, -5, -5)), max(
				coordinate(5, 5, 5)) {

	// think about parameteres for the individual pattern... give as parameters with default value

	switch (pattern) {
	case 1: { // poisson pattern
		std::cout
				<< "Generate Poisson pattern with 1000 points in 10^3 cubicle..."
				<< std::endl;

		name = "Poisson";

		// seed for random numbers generator
		srand(time(NULL));

		for (unsigned n = 0; n < 100; n++) {
			list.push_back(
					new coordinate(coordinate(3) * 10 - getExtend() / 2));
		}

		break;
	}
	case 2: { // diamond point pattern ((c) Dirk Ropers)
		std::cout
				<< "Generate Diamond point pattern with 1000 points in 10^3 cubicle..."
				<< std::endl;

		name = "Diamond";

		// TODO: generate unit cell and then extend to desired size
		// TODO: this must be the box size? Gives point density of one...
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
				<< "Generate Test point pattern with explicitly given points..."
				<< std::endl;

		name = "Test";

		// TODO: Add points per gui?

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

		min = coordinate(-2.5, -2.5, -2.5);
		max = coordinate(2.5, 2.5, 2.5);

		break;
	}
	}

	std::cout << name << " point pattern generated." << std::endl;
}

Point3 pointlist::getPoint(const Periodic_point& pt, const PDT& PD3d) const {
	CGAL::Periodic_3_offset_3 offset = pt.second;

	for (unsigned int i = 0; i < 3; i++) {

		// shift every offset larger than 2 back
		// 				  n = 0 1 2  3 4 5
		// (((n+1)mod 3)-1) = 0 1 -1 0 1 -1

		offset[i] = (((offset[i] + 1) % 3) - 1);
	}

	return PD3d.construct_point(pt.first, offset);
}

int pointlist::pointsInside(const std::vector<coordinate>& points,
		const coordinate& mid, const double r) const {
	if (r == 0) {
		return 0;
	}
	int ctr = 0;

	// nodesiterations
	for (std::vector<coordinate>::const_iterator point = points.begin();
			point != points.end(); ++point) {
		// bounding box of the sphere
		if ((fabs((*point)[0] - mid[0]) < r) && (fabs((*point)[1] - mid[1]) < r)
				&& (fabs((*point)[2] - mid[2]) < r)) {
			// TODO: more elegant but slower by a factor of 3. why?
			//if ((*point).insideAABB(mid, r)) {
			// sphere itself
			if (((*point) - mid).lengthSqr() < r * r) {
				ctr++;
			}
		}
	}
	return ctr;
}

void pointlist::add(coordinate pos) {
	add(pos.x(), pos.y(), pos.z());
}

void pointlist::add(double x, double y, double z) {
	coordinate* newPos = new coordinate(x, y, z);

	// TODO: same as getShifter?! maybe put function to shift point to domain in coordinate...
	// ensure that point lies in between boundaries if periodic
	if (periodic) {
		coordinate extend = getExtend();
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
/*
 void pointlist::setExtend(coordinate ext) {
 extend = ext;
 }
 */
coordinate pointlist::getExtend() const {
	return max - min;
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
	// correct domain size
	min *= a;
	max *= a;

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

	// initialise nodelist
	nodelist* nlist = new nodelist(periodic, name, this);
	nlist->setMins(min);
	nlist->setMaxs(max);

	// **************** debugging *************
	/*
	 // initialize geomview for visulaization
	 CGAL::Geomview_stream gv(
	 CGAL::Bbox_3(min.x(), min.y(), min.z(), max.x(), max.y(), max.z()));
	 gv.set_line_width(4);
	 // gv.set_trace(true);
	 gv.set_bg_color(CGAL::Color(0, 200, 200));
	 gv.clear();
	 */
	// **************** debugging *************
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

			// **************** debugging *************
			/*
			 // verification
			 if (list.size() < 100) {
			 gv
			 << Sphere3(Point3((*c)->x(), (*c)->y(), (*c)->z()),
			 (double) 0.005);
			 }
			 */
			// **************** debugging *************
			++show_progress;
		}

		// iteration over unique tetrahedra that intersects the original domain (centre may lie in between boundaries even if corners of tetrahedron lies outside)
		// initialisation of objects
		PDT::Periodic_tetrahedron pt;
		Tetrahedron tet;
		PDT::Cell_handle cell;
		PDT::Cell_handle neighbourCell;

		Point3 nodePoint;
		coordinate nodeCentre;
		node* n;

		Point3 neighbourPoint;
		coordinate neighbourCentre;
		node* neighbour;

		std::cout << "is onesheeted:" << PD3d.is_triangulation_in_1_sheet()
				<< std::endl;

		// estimate duration
		std::cout
				<< "Calculate centroids (nodes) and links to attached tetrahedra (to neighbours)..."
				<< std::endl;
		show_progress.restart(
				std::distance(
						PD3d.periodic_tetrahedra_begin(
								PDT::UNIQUE_COVER_DOMAIN),
						PD3d.periodic_tetrahedra_end(
								PDT::UNIQUE_COVER_DOMAIN)));

		// TODO: speed up the process by iterating through cells and skip the locate part?
		for (PDT::Periodic_tetrahedron_iterator ptit =
				PD3d.periodic_tetrahedra_begin(PDT::UNIQUE_COVER_DOMAIN);
				ptit != PD3d.periodic_tetrahedra_end(PDT::UNIQUE_COVER_DOMAIN);
				++ptit) {
			// get tetrahedron
			pt = *ptit;

			// Convert the current Periodic_tetrahedron in T^3 to a Tetrahedron in R^3
			tet = PD3d.construct_tetrahedron(pt);

			// calculate centre
			nodePoint = CGAL::centroid(tet.vertex(0), tet.vertex(1),
					tet.vertex(2), tet.vertex(3));
			nodeCentre = coordinate(nodePoint.x(), nodePoint.y(),
					nodePoint.z());

			// to get the progressbar right before skipping a cell...
			++show_progress;

			// don't evaluate tetrahedron if its centroid is not in the original volume
			if (!nodeCentre.insideAABB(min, max)) {
				continue;
			}

			// add node to nodelist
			n = nlist->add(nodeCentre);

			// get cell for the current nodepoint. Most timeconsuming part?
			cell = PD3d.locate(nodePoint); // has an offset
			//cell = ptit.get_cell(); // does not have an offset

			// **************** debugging *************
			/*
			 //std::cout << "cell offset 0 " << cell->offset(0) << std::endl;
			 //std::cout << "cell offset 1 " << cell->offset(1) << std::endl;
			 //std::cout << "cell offset 2 " << cell->offset(2) << std::endl;
			 //std::cout << "cell offset 3 " << cell->offset(3) << std::endl;


			 std::cout << "Iterated tetrahedron + centre" << std::endl;
			 gv << CGAL::RED;
			 gv << Sphere3(tet.vertex(0), 0.01);
			 gv << Sphere3(tet.vertex(1), 0.01);
			 gv << Sphere3(tet.vertex(2), 0.01);
			 gv << Sphere3(tet.vertex(3), 0.01);
			 gv << tet;
			 gv << CGAL::GREEN;
			 gv << Sphere3(nodePoint, 0.2);

			 std::cout << "Continue with Enter." << std::endl;
			 std::cin.get();
			 */
			// **************** debugging *************
			// neighbouriteration
			for (int i = 0; i < 4; i++) {
				// get current neighbouring cell
				neighbourCell = cell->neighbor(i);

				// construct neighbouring tetrahedron midpoint (by applying the periodicity as in the periodicpointpatterntool)
				neighbourPoint = CGAL::centroid(
						getPoint(PD3d.periodic_point(neighbourCell->vertex(0)),
								PD3d),
						getPoint(PD3d.periodic_point(neighbourCell->vertex(1)),
								PD3d),
						getPoint(PD3d.periodic_point(neighbourCell->vertex(2)),
								PD3d),
						getPoint(PD3d.periodic_point(neighbourCell->vertex(3)),
								PD3d));
				neighbourCentre = coordinate(neighbourPoint.x(),
						neighbourPoint.y(), neighbourPoint.z());

				// add the neighbouring node
				neighbour = nlist->add(neighbourCentre);

				// TODO: put function addPair(coordinate p1, coordinate p2) in nodelist?
				// set neighbourhood
				n->addNeighbour(neighbour, nlist->getShifter(nodeCentre),
						nlist->getShifter(neighbourCentre));
				neighbour->addNeighbour(n, nlist->getShifter(neighbourCentre),
						nlist->getShifter(nodeCentre));

				// **************** debugging *************
				/*
				 std::cout << "neighbour " << i << "tetrahedron + centre from cell" << std::endl;
				 //std::cout << "neighbourCell offset 0 " << neighbourCell->offset(0) << std::endl;
				 //std::cout << "neighbourCell offset 1 " << neighbourCell->offset(1) << std::endl;
				 //std::cout << "neighbourCell offset 2 " << neighbourCell->offset(2) << std::endl;
				 //std::cout << "neighbourCell offset 3 " << neighbourCell->offset(3) << std::endl;
				 gv << CGAL::YELLOW;
				 gv << Sphere3(getPoint(PD3d.periodic_point(neighbourCell->vertex(0)), PD3d),0.02);
				 gv << Sphere3(getPoint(PD3d.periodic_point(neighbourCell->vertex(1)), PD3d),0.02);
				 gv << Sphere3(getPoint(PD3d.periodic_point(neighbourCell->vertex(2)), PD3d),0.02);
				 gv << Sphere3(getPoint(PD3d.periodic_point(neighbourCell->vertex(3)), PD3d),0.02);

				 Tetrahedron neighTet = PD3d.construct_tetrahedron(PD3d.construct_periodic_3_tetrahedron(getPoint(PD3d.periodic_point(neighbourCell->vertex(0)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(1)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(2)), PD3d), getPoint(PD3d.periodic_point(neighbourCell->vertex(3)), PD3d)));
				 gv << neighTet;

				 gv << CGAL::BLUE;
				 gv << Sphere3(neighbourPoint,0.2);

				 std::cout << "Continue with Enter." << std::endl;
				 std::cin.get();
				 */
				// **************** debugging *************
			}

			// **************** debugging *************
			//gv.clear();
			// **************** debugging *************
		}
	} else {
		// use http://doc.cgal.org/latest/Triangulation_3/index.html

		// calculate delaunay triangulation
		DT D3d;
		std::vector<Point3> points(this->size());
		for (std::vector<coordinate*>::iterator c = list.begin();
				c != list.end(); ++c) {
			D3d.insert(Point3((*c)->x(), (*c)->y(), (*c)->z()));
			// **************** debugging *************
			/*
			 if (list.size() < 100) {
			 gv
			 << Sphere3(Point3((*c)->x(), (*c)->y(), (*c)->z()), (double) 0.1);
			 }
			 */
			// **************** debugging *************
		}
		// **************** debugging *************
		/*
		 gv.set_wired(true);
		 gv << D3d;
		 gv << CGAL::GREEN;
		 D3d.draw_dual(gv);
		 gv.set_wired(false);
		 */
		// **************** debugging *************
		boost::progress_display show_progress(D3d.number_of_cells());

		//int i = 0;
		// Iterate over every tetrahedra ("cell" internally)
		for (DT::Finite_cells_iterator cell = D3d.finite_cells_begin();
				cell != D3d.finite_cells_end(); ++cell) {

			// calculate midpoint and add to nodelist
			Point3 centre = CGAL::centroid((*cell).vertex(0)->point(),
					(*cell).vertex(1)->point(), (*cell).vertex(2)->point(),
					(*cell).vertex(3)->point());

			node* n = nlist->add(centre.x(), centre.y(), centre.z());
			// **************** debugging *************
			/*
			 gv << CGAL::RED;
			 gv << Sphere3(centre, (double) 0.005);
			 */
			// **************** debugging *************
			// iterate of all four neighbours, calculate their centroid and set the neighbourhood. TODO: points are being added fourfold...
			for (int neighIdx = 0; neighIdx < 4; neighIdx++) {
				// exclude infinite cells
				if (!D3d.is_infinite((*cell).neighbor(neighIdx))) {

					Point3 neighbour = CGAL::centroid(
							(*cell).neighbor(neighIdx)->vertex(0)->point(),
							(*cell).neighbor(neighIdx)->vertex(1)->point(),
							(*cell).neighbor(neighIdx)->vertex(2)->point(),
							(*cell).neighbor(neighIdx)->vertex(3)->point());
					// **************** debugging *************
					/*
					 gv << CGAL::BLUE;
					 gv << Sphere3(neighbour, (double) 0.005);
					 */
					// **************** debugging *************
					node* neigh = nlist->add(neighbour.x(), neighbour.y(),
							neighbour.z());

					// set the neighbourhood
					n->addNeighbour(neigh);
					neigh->addNeighbour(n);
				}
			}

			++show_progress;
		}
		// set edgenodes
		double characteristicLength = stats(nlist->lengthDistribution())[1] * 2;
		std::cout << "Nodes farther away of the edge than "
				<< characteristicLength << " are declared to be edgenodes."
				<< std::endl;
		nlist->setEdgenodes(characteristicLength);
	}
	t = clock() - t;
	std::cout << "Took " << (((float) t) / CLOCKS_PER_SEC)
			<< "s to generate the neighbours." << std::endl;
	// **************** debugging *************

	/*
	 std::cout << "Continue with Enter." << std::endl;
	 std::cin.get();*/
	// **************** debugging *************
	return nlist;
}

std::vector<std::vector<double> > pointlist::hyperuniformity(unsigned int nr,
		unsigned int n) const {
	// radiusincrement + maximal radius
	double rMax = periodic ? getExtend().min() : getExtend().min() / 3, dr =
			rMax / nr;

	// datastructure to save the number of points in each sphere
	std::vector<std::vector<double> > data;
	data.resize(nr);
	for (unsigned int i = 0; i < nr; ++i) {
		data[i].resize(n);
		// do not calculate radius=0
		data[0][i] = 0;
	}

	clock_t t;
	t = clock();

	// seed for generation of random numbers
	srand(time(NULL));

	// midpoint of sphere
	coordinate mid;

	if (periodic) {
		// precalculate extended pattern
		std::vector<coordinate> shifters = getExtend().getShifters();
		std::vector<coordinate> extendedPattern;
		for (std::vector<coordinate>::iterator shifter = shifters.begin();
				shifter != shifters.end(); ++shifter) {
			for (std::vector<coordinate*>::const_iterator n = list.begin();
					n != list.end(); ++n) {
				extendedPattern.push_back(*(*n) + (*shifter));
			}
		}

		// iteration over n spheres
		boost::progress_display show_progress(n);
		for (unsigned int i = 0; i < n; i++) {
			// Choose the centre of the sphere such that it is somewhere within the whole pattern.
			mid = (coordinate(3) - 0.5) * getExtend();

			// iteration over radius: count points within sphere
			for (unsigned int j = 1; j < nr; j++) {
				data[j][i] = pointsInside(extendedPattern, mid, j * dr);
			}

			++show_progress;
		}
	} else {
		// get patterns coordinates
		std::vector<coordinate> pattern;
		for (std::vector<coordinate*>::const_iterator n = list.begin();
				n != list.end(); ++n) {
			pattern.push_back(**n);
		}

		double r; //, rSqr;
		// iteration over radius
		boost::progress_display show_progress(nr);
		for (unsigned int j = 0; j < nr; j++) {
			r = j * dr;
			//rSqr = r * r;
			// iteration over n spheres
			for (unsigned int i = 0; i < n; i++) {
				// Choose the center of the sphere such that it is somewhere within the whole pattern.
				mid = (coordinate(3) - 0.5) * (getExtend() - 2 * r);

				data[j][i] = pointsInside(pattern, mid, r);
			}
			++show_progress;
		}
	}
	t = clock() - t;
	std::cout << "It took " << (((float) t) / CLOCKS_PER_SEC)
			<< "s to calculate the hyperuniformity of pattern " << name << "."
			<< std::endl;

	//cout << "Radius" << char(9) << "Volumen" << char(9) << "Erwartungswert" << char(9) << "Verhältnis" << endl;
	// Calculate the expected value (of numbers of points within sphere) of each radius. Should roughly be equal to the numberdensity times volume of the sphere.
	std::vector<double> expectedValue;
	expectedValue.resize(nr);
	// iteration over radius
	for (unsigned int j = 0; j < nr; j++) {
		// iteration over sphere
		for (unsigned int i = 0; i < n; i++) {
			expectedValue[j] += data[j][i];
		}
		expectedValue[j] = expectedValue[j] / n;
		//cout << dr*j << char(9) << 4/3 * M_PI * pow(dr*j, 3) << char(9) << expectedValue[j] << char(9) << expectedValue[j]/(4/3 * M_PI * pow(dr*j, 3)) << endl;
	}

	// Calculate variance for each radius. variance[i][j], (i,j) = (rows,colums) = (radius, variance(radius))
	std::vector<std::vector<double> > variance;
	variance.resize(nr);
	// iteration over radius
	for (unsigned int j = 0; j < nr; j++) {
		variance[j].resize(2);
		variance[j][0] = j * dr;

		// iteration over sphere
		for (unsigned int i = 0; i < n; i++) {
			variance[j][1] += pow((data[j][i] - expectedValue[j]), 2);
		}
		variance[j][1] = variance[j][1] / n;
	}

	return variance;
}

std::vector<std::vector<double> > pointlist::structurefactor(
		const unsigned int nq) const {
	// TODO: generate lists of wave vectors q, randomly distributed on spheres of different radii
	std::vector<std::vector<std::vector<double> > > qs;

	// TODO: calculate sum exp(iqr) for each q vector and average value for qs of same length:

	std::vector<std::vector<double> > sq;
	sq.resize(nq);
	// iteration over radius
	for (unsigned int q = 0; q < nq; q++) {
		sq[q].resize(2);
	}

	// test data for plotting

	sq[0][0] = 0; // q
	sq[0][1] = 0; // S(q)
	sq[1][0] = 1; // q
	sq[1][1] = 0.5; // S(q)
	sq[2][0] = 2; // q
	sq[2][1] = 0.5; // S(q)
	sq[3][0] = 3; // q
	sq[3][1] = 1; // S(q)
	sq[4][0] = 4; // q
	sq[4][1] = 5; // S(q)

	std::cout << "not yet implemented, return dummy function" << std::endl;

	return sq;
}
