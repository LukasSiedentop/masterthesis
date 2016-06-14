/*
 * nodelist.cpp
 *
 *  Created on: 24.07.2015
 *      Author: Lukas Siedentop
 */

#include "nodelist.hpp"

//using namespace std;

nodelist::nodelist() :
		periodic(0), name(""), min(coordinate()), max(coordinate()) {
}

nodelist::nodelist(bool periodicity, std::string n) :
		periodic(periodicity), name(n), min(coordinate()), max(coordinate()) {
}

nodelist::nodelist(int pattern, bool periodicity) :
		periodic(periodicity), min(coordinate(-5, -5, -5)), max(
				coordinate(5, 5, 5)) {
	std::cout
			<< "Warning! Decorated pattern is not created via underlying point pattern but directly."
			<< std::endl;
	switch (pattern) {
	case 1: { // poisson pattern
		std::cout
				<< "Generate Poisson pattern with 1000 points in 10^3 cubicle..."
				<< std::endl;

		name = "Poisson";

		// seed for random numbers generator
		srand(time(NULL));
		for (unsigned n = 0; n < 1000; n++) {
			node* tmp = new node(this,
					(coordinate(3) * 10) - coordinate(5, 5, 5));
			tmp->setEdgenode(1);
			list.push_back(tmp);
		}

		break;
	}
	case 2: { // diamond point pattern ((c) Dirk)
		std::cout
				<< "Generate Diamond point pattern with 1000 points in 10^3 cubicle..."
				<< std::endl;

		name = "Diamond";

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
								new node(this,
										(coordinate(0.25, 0.25, 0.25)
												+ coordinate(x, y, z))));
						// second point: (0.75,0.75,0.75)
						list.push_back(
								new node(this,
										(coordinate(0.75, 0.75, 0.75)
												+ coordinate(x, y, z))));
					}
					x++;
				}
				y++;
			}
			z++;
		}

		// set edgenodes
		for (std::vector<node*>::iterator n = list.begin(); n != list.end();
				++n) {
			(*n)->setEdgenode(1);
		}

		break;
	}

	case 3: { // test pattern (tetrahedron
		std::cout
				<< "Generate Test point pattern with 4 points in 10^3 cubicle..."
				<< std::endl;

		name = "Test";
		/*
		 // Tetrahedron
		 list.push_back(new node(this, (coordinate(-1, 0, -1))));
		 list.push_back(new node(this, (coordinate(1, 0, -1))));
		 list.push_back(new node(this, (coordinate(0, -1, 1))));
		 list.push_back(new node(this, (coordinate(0, 1, 1))));
		 */

		list.push_back(new node(this, (coordinate(0.1, 0, 0))));
		list.push_back(new node(this, (coordinate(-1.2, -1, -1))));
		list.push_back(new node(this, (coordinate(-1.3, -1, 1))));
		list.push_back(new node(this, (coordinate(-1.4, 1, -1))));
		list.push_back(new node(this, (coordinate(-1.5, 1, 1))));
		list.push_back(new node(this, (coordinate(1.1, -1, -1))));
		list.push_back(new node(this, (coordinate(1.2, -1, 1))));
		list.push_back(new node(this, (coordinate(1.3, 1, -1))));
		list.push_back(new node(this, (coordinate(1.4, 1, 1))));

		min = coordinate(-1.5, -1.5, -1.5);
		max = coordinate(1.5, 1.5, 1.5);

		break;
	}
	}
	setNeighbours();

	std::cout << name << " point pattern generated. Statistics:" << std::endl
			<< listStats() << std::endl;
}

bool nodelist::isPeriodic() {
	return periodic;
}

void nodelist::setNeighbours(unsigned int valency) {
	// TODO: may be computationally optimized: outer loop probably too large
	// Method should give similar results as designprotocol for crystals

	std::cout << "Generating the neighbours for the pattern " << name
			<< " in a very unsophisticated way (the closest four without four neighbours)."
			<< std::endl;

	// meassure the time
	clock_t t;
	t = clock();

	boost::progress_display show_progress(valency * list.size());
	// get four next points
	for (unsigned int neighs = 0; neighs < valency; neighs++) {
		// nodesiteration
		for (std::vector<node*>::iterator n = list.begin(); n != list.end();
				++n) {
			node* neigh = NULL;
			// if selected node has less than 4 (valency) neighbours
			if ((*n)->countNeighbours() < valency) {

				// find closest neighbour with less than four neighbours
				double distSqr = std::numeric_limits<double>::infinity();
				for (std::vector<node*>::iterator nextn = list.begin();
						nextn != list.end(); ++nextn) {

					// don't count if its already a neighbour, is the node itself or valency is already satisfied
					if ((*n)->isNeighbour((*nextn)) || n == nextn
							|| (*nextn)->countNeighbours() >= valency) {
						continue;
					}

					// connecting vector of node and going to be neighbour
					coordinate diff;

					// compare the difference vectors with an added shifter for the periodic case
					if (periodic) {
						std::vector<coordinate> shifters = getShifters();
						for (std::vector<coordinate>::iterator shifter =
								shifters.begin(); shifter != shifters.end();
								++shifter) {

							// calculate connecting vector
							diff = (((*nextn)->getPosition()) + (*shifter))
									- ((*n)->getPosition());

							// check if the connection is shorter than the previously saved one
							if (diff.lengthSqr() < distSqr) {
								distSqr = diff.lengthSqr();
								neigh = (*nextn);
							}
						}
					} else {
						// calculate connecting vector
						diff = ((*nextn)->getPosition())
								- ((*n)->getPosition());

						// check if the connection is shorter than the previously saved one
						if (diff.lengthSqr() < distSqr) {
							distSqr = diff.lengthSqr();
							neigh = (*nextn);
						}
					}
				}
			}

			if (neigh) {
				(*n)->addNeighbour(neigh);
				neigh->addNeighbour(*n);
			}

			++show_progress;
		}
	}

	t = clock() - t;
	std::cout << "Took " << (((float) t) / CLOCKS_PER_SEC)
			<< "s to generate the neighbours." << std::endl;
}

void nodelist::setDensity(double density) {
	scaleList(pow((getDensity() / density), 1. / 3.));
}

double nodelist::getDensity() {
	return list.size() / getVolume();
}

double nodelist::getVolume() {
	double volume = 1;

	coordinate lengths = getLengths();

	for (unsigned int i = 0; i < lengths.dimensions(); i++) {
		volume *= lengths[i];
	}

	return volume;
}

int nodelist::countEdgenodes() {
	if (periodic) {
		return 0;
	}
	int ctr = 0;
	for (std::vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		if ((*n)->isEdgenode()) {
			ctr++;
		}
	}
	return ctr;
}

double nodelist::r(double theta, double w, double h) {
	return (w * h)
			/ (2
					* sqrt(
							w * w * sin(theta) * sin(theta)
									+ h * h * cos(theta) * cos(theta)));
}

int nodelist::pointsInside(const std::vector<coordinate>& points,
		const coordinate& mid, const double r) const {
	if (r==0){
		return 0;
	}
	int ctr = 0;

	// nodesiterations
	for (std::vector<coordinate>::const_iterator point = points.begin();
			point != points.end(); ++point) {


		// bounding box of the sphere

		if ((fabs((*point)[0] - mid[0]) < r) && (fabs((*point)[1] - mid[1]) < r)
				&& (fabs((*point)[2] - mid[2]) < r)) {
		/*
		std::cout << (mid-r).toString() << std::endl;
		std::cout << (mid+r).toString() << std::endl;
		if ((*point).insideAABB(mid-r, mid+r)) {*/
			// sphere itself
			if (((*point) - mid).lengthSqr() < r * r) {
				ctr++;
			}
		}
	}
	return ctr;
}
/*
 nodelist nodelist::extendList(int nx, int ny, int nz) {
 double lx = getLengths()[0], ly = getLengths()[1], lz = getLengths()[2];

 nodelist* extList = new nodelist(false, name + " extended x " + nx + ", y " + ny + ", z " + nz);

 // go through all combinations
 for (int ix = 0; ix < nx; ix++) {
 for (int iy = 0; iy < ny; iy++) {
 for (int iz = 0; iz < nz; iz++) {
 // add the shifting vector to the list and merge it
 coordinate(ix * lx, iy * ly, iz * lz);
 extList->merge(this->shiftList());
 }
 }
 }

 nodelist* extList = new nodelist(this->shiftList(shifters[i]));



 return this;
 }
 */
void nodelist::setMins(coordinate mins) {
	min = mins;
}

void nodelist::setMaxs(coordinate maxs) {
	max = maxs;
}

void nodelist::shiftList(coordinate shifter) {
	// correct bounding box
	min += shifter;
	max += shifter;

	// shift entries
	for (std::vector<node*>::iterator it = list.begin(); it != list.end();
			++it) {
		(*it)->shift(shifter);
	}
}

void nodelist::scaleList(double a) {
	// correct bounding box
	min *= a;
	max *= a;

	// scale entries
	for (std::vector<node*>::iterator it = list.begin(); it != list.end();
			++it) {
		(*it)->scale(a);
	}
}

void nodelist::scaleListAnisotropic(double ax, double ay, double az) {
	// correct bounding box
	min[0] *= ax;
	min[1] *= ay;
	min[2] *= az;
	max[0] *= ax;
	max[1] *= ay;
	max[2] *= az;

	// scale entries
	for (std::vector<node*>::iterator it = list.begin(); it != list.end();
			++it) {
		(*it)->scaleAnisotropic(ax, ay, az);
	}
}

/* unneccesary
 bool nodelist::isPeriodic() const {
 return periodic;
 }
 */

std::string nodelist::getName() {
	return name;
}

node* nodelist::add(coordinate pos) {
	return add(pos.x(), pos.y(), pos.z());
}

node* nodelist::add(double x, double y, double z) {
	coordinate newPos = coordinate(x, y, z);

	// ensure that point lies in between boundaries if periodic
	if (periodic) {
		std::vector<coordinate> shifters = getShifters();
		for (std::vector<coordinate>::iterator shifter = shifters.begin();
				shifter != shifters.end(); ++shifter) {
			// as bounding box: choose the right shifter where all components are in t
			/*
			if ((fabs((*shifter).x() + (x + getMid().x()))
					< getLengths().x() / 2)
					&& (fabs((*shifter).y() + (y + getMid().y()))
							< getLengths().y() / 2)
					&& (fabs((*shifter).z() + (z + getMid().z()))
							< getLengths().z() / 2)) {*/
			if ((newPos+(*shifter)).insideAABB(min, max)) {


				newPos = (*shifter) + coordinate(x, y, z);
				break;
			}
		}
	}

	// check if point at given position exist
	for (std::vector<node*>::iterator it = list.begin(); it != list.end();
			++it) {
		if ((*it)->getPosition() == newPos) {
			return *it;
		}
	}

	// if not, build and add new node
	node* n = new node(this, newPos);
	list.push_back(n);
	return n;
}

 coordinate nodelist::getMins() {
 return min;
 }

 coordinate nodelist::getMaxs() {
 return max;
 }

coordinate nodelist::getLengths() {
	return max - min;
}

coordinate nodelist::getMid() {
	return (max + min) / 2;
}

double nodelist::getMaxFeatureSize() {
	return getLengths().length() / 2;
}

void nodelist::display() {
	std::cout << "Number of elements: " << list.size() << std::endl;

	int i = 1;

	// nodesiteration
	for (std::vector<node*>::iterator it = list.begin(); it != list.end();
			++it) {
		std::cout << char(9) << i << "-th element: " << ((*it)->getPosition())
				<< " with " << (*it)->countNeighbours() << " neighbours at: ";
		// neighbouriteration
		for (std::vector<class node *>::iterator neighIt =
				(*it)->getNeighbours()->begin();
				neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
			std::cout << ((*neighIt)->getPosition()) << " ";
		}
		std::cout << std::endl;

		i++;
	}
}

coordinate nodelist::getShifter(coordinate point) {
	if (!periodic) {
		std::cout
				<< "WARNING: shifter vector requested whilst the pattern is not periodic!";
		return coordinate(0, 0, 0);
	}

	// ensure that point lies in between boundaries if periodic
	std::vector<coordinate> shifters = getShifters();

	unsigned int nBoxes = 100;
	for (unsigned int n = 1; n < nBoxes; n++) {
		if (n > 1) {
			std::cout
					<< "Note: No Shifter vector found; requested point may be "
					<< n << " boxes away..." << std::endl;
		}

		for (std::vector<coordinate>::iterator shifter = shifters.begin();
				shifter != shifters.end(); ++shifter) {
			coordinate largeShifter = (*shifter) * n;
			// as bounding box: choose the right shifter where all components are in t
			if ((fabs(largeShifter.x() + (point.x() + getMid().x()))
					< getLengths().x() / 2)
					&& (fabs(largeShifter.y() + (point.y() + getMid().y()))
							< getLengths().y() / 2)
					&& (fabs(largeShifter.z() + (point.z() + getMid().z()))
							< getLengths().z() / 2)) {
				return largeShifter;
			}
		}
	}
	std::cout << "Note: no shifting vector was found for point "
			<< point.toString() << " into the box " << min.toString() << ", "
			<< max.toString() << " despite looking for shifters over " << nBoxes
			<< " boxes. Return zero Shifter." << std::endl;

	return coordinate(0, 0, 0);
}

std::vector<coordinate> nodelist::getShifters() {
	return getLengths().getShifters();
}

std::vector<coordinate> nodelist::getShifted(coordinate mid,
		double halfExtend) {
	std::vector<coordinate> shifted;

	// boxlength
	double lx = getLengths()[0], ly = getLengths()[1], lz = getLengths()[2];

	// go through all combinations, but only if one edge of the given box is not within the pattern...
	for (int ix = -1; ix < 2; ix++) {
		if ((fabs(mid[0] + (ix * halfExtend)) >= fabs(ix * (lx / 2.0)))) {
			for (int iy = -1; iy < 2; iy++) {
				if ((fabs(mid[1] + (iy * halfExtend)) >= fabs(iy * (ly / 2.0)))) {
					for (int iz = -1; iz < 2; iz++) {
						if ((fabs(mid[2] + (iz * halfExtend))
								>= fabs(iz * (lz / 2.0)))) {
							// ... add the shifting vector (Shift sphere contrary to pattern!)
							shifted.push_back(
									mid
											- coordinate(ix * lx, iy * ly,
													iz * lz));
						}
					}
				}
			}
		}
	}
	return shifted;
}

std::string nodelist::listStats(const std::string commentDelimeter) {
	std::stringstream stream;
	stream << commentDelimeter << (periodic ? "Periodic" : "Non-periodic")
			<< std::endl;
	if (!periodic) {
		stream << commentDelimeter << "Number of edgenodes: "
				<< countEdgenodes() << std::endl;
	}

	stream << commentDelimeter << "Number of nodes: " << list.size()
			<< std::endl;
	stream << commentDelimeter << "Bounding box: " << min << ", " << max
			<< std::endl;

	stream << commentDelimeter << "Characteristic length: "
			<< stats(lengthDistribution())[1] << std::endl;

	stream << commentDelimeter << "Mid of box: " << getMid() << std::endl;
	stream << commentDelimeter << "Volume: " << getVolume() << std::endl;
	stream << commentDelimeter << "Density of points: " << getDensity()
			<< std::endl;

	return stream.str();
}

double nodelist::normalize() {
	// stats(lengthDistribution())[1]

	// only for agapornis! TODO: function to set scaling for each direction...
	/*
	 double resX = 9.067; // nm/px
	 double resY = 9.067; // nm/px
	 double resZ = 25; //nm/px

	 scaleListAnisotropic(resX, resY, resZ);
	 */

	// midpoint
	shiftList(getMid() * -1);

	double factor = pow(getDensity(), 1. / 3.);

	// set density to 1
	setDensity(1);

	return factor;
}

std::vector<double> nodelist::neighbourDistribution() {
	std::vector<double> data;

	int counter;

	// nodesiteration
	for (std::vector<node*>::iterator it = list.begin(); it != list.end();
			++it) {
		if (!(*it)->isEdgenode()) {
			counter = 0;
			// neighboursiteration
			for (std::vector<node*>::iterator neighIt =
					(*it)->getNeighbours()->begin();
					neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
				counter++;
			}

			data.push_back(counter);
		}
	}

	return data;
}

std::vector<double> nodelist::lengthDistribution() {
	std::vector<double> data;

	// nodesiteration
	for (std::vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		// exclude distances between edgendoes
		if (!(*n)->isEdgenode()) {
			// neighboursiteration
			for (std::vector<node*>::iterator nn =
					(*n)->getNeighbours()->begin();
					nn != (*n)->getNeighbours()->end(); ++nn) {
				if (periodic) {
					data.push_back((*n)->euklidianPeriodic((*nn)));
				} else {
					data.push_back((*n)->euklidian((*nn)));
					// if the neighbour is an edgenode, put the distance twice, as later every second distance is taken
					if ((*nn)->isEdgenode()) {
						data.push_back((*n)->euklidian((*nn)));
					}
				}
			}
		}
	}

	// sort data and keep only every second
	data = mergeSort(data);
	std::vector<double> halfData;
	for (unsigned int i = 0; i < data.size(); i += 2) {
		halfData.push_back(data[i]);
	}

	return halfData;
}

std::vector<double> nodelist::angleDistribution() {
	std::vector<double> data;

	// nodesiteration
	for (std::vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		if (!(*nodeIter)->isEdgenode()) {
			// neighboursiteration
			for (std::vector<node*>::iterator neighIter1 =
					(*nodeIter)->getNeighbours()->begin();
					neighIter1 != (*nodeIter)->getNeighbours()->end();
					++neighIter1) {

				// neighbours of neighboursiteration
				for (std::vector<node*>::iterator neighIter2 =
						(*neighIter1)->getNeighbours()->begin();
						neighIter2 != (*neighIter1)->getNeighbours()->end();
						++neighIter2) {
					// don't calculate angles with itself
					if (!(*nodeIter)->equals(*neighIter2)) {
						if (periodic) {
							data.push_back(
									((180.0) / M_PI)
											* (*neighIter1)->anglePeriodic(
													(*nodeIter),
													(*neighIter2)));
						} else {

							data.push_back(
									((180.0) / M_PI)
											* (*neighIter1)->angle((*nodeIter),
													(*neighIter2)));
						}
					}
				}
			}
		}
	}

	// sort data and keep only every second
	data = mergeSort(data);
	std::vector<double> halfData;
	for (unsigned int i = 0; i < data.size(); i += 2) {
		halfData.push_back(data[i]);
	}

	return halfData;
}

std::vector<std::vector<double> > nodelist::hyperuniformity(unsigned int nr,
		unsigned int n) {
	// radiusincrement + maximal radius
	double rMax = periodic ? getLengths().min() : getLengths().min() / 3, dr =
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
		std::vector<coordinate> shifters = getShifters();
		std::vector<coordinate> extendedPattern;
		for (std::vector<coordinate>::iterator shifter = shifters.begin();
				shifter != shifters.end(); ++shifter) {
			for (std::vector<node*>::iterator n = list.begin(); n != list.end();
					++n) {
				extendedPattern.push_back((*n)->getPosition() + *shifter);
			}
		}

		// iteration over n spheres
		boost::progress_display show_progress(n);
		for (unsigned int i = 0; i < n; i++) {
			// Choose the centre of the sphere such that it is somewhere within the whole pattern.
			mid = (coordinate(3) - 0.5) * getLengths();

			// iteration over radius: count points within sphere
			for (unsigned int j = 1; j < nr; j++) {
				data[j][i] = pointsInside(extendedPattern, mid, j * dr);
			}

			++show_progress;
		}
	} else {
		// get patterns coordinates
		std::vector<coordinate> pattern;
		for (std::vector<node*>::iterator n = list.begin(); n != list.end();
				++n) {
			pattern.push_back((*n)->getPosition());
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
				mid = (coordinate(3) - 0.5) * (getLengths() - 2 * r);

				data[j][i] = pointsInside(pattern, mid, r);
			}
			++show_progress;
		}
	}
	t = clock() - t;
	std::cout << "Took " << (((float) t) / CLOCKS_PER_SEC) << "s" << std::endl;

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

void nodelist::writeCoordinates() {
	std::stringstream fileNameStream;
	fileNameStream << "./data/points_" << name << ".csv";
	std::string outfileName = fileNameStream.str();

	std::vector<std::string> data;
	/*
	 // iterate over extended pattern
	 std::vector<coordinate> shifters = getShifters();

	 for (vector<coordinate>::iterator shifter = shifters.begin();
	 shifter != shifters.end(); ++shifter) {
	 for (vector<node*>::iterator n = list.begin(); n != list.end();
	 ++n) {
	 stringstream stream;
	 stream << ((*n)->getPosition() + *shifter).toString("", "\t", "") << endl;
	 data.push_back(stream.str());
	 }
	 }
	 */

	// iterate over pattern
	for (std::vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		std::stringstream stream;
		stream << ((*n)->getPosition()).toString("", "\t", "", 10) << std::endl;
		data.push_back(stream.str());
	}

	// write to outfile
	std::ofstream outfile;
	outfile.open(outfileName.c_str());
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	std::cout << "Pattern written to file '" << outfileName << "'. Enjoy!"
			<< std::endl;
}

void nodelist::writeGWL() {
	// TODO: work on list that is repeated to the wanted size
	int x = 2;
	int y = 2;
	int z = 1;
	std::vector<class node*> extendedList((x + y + z) * size());

	// sort by metric
	// TODO: make sortMetric chooseable (cone, ...)
	sort(extendedList.begin(), extendedList.end(), node::zmetric);

	// get loops
	std::vector<coordinate> loop;
	for (std::vector<node*>::iterator nodeIter = extendedList.begin();
			nodeIter != extendedList.end(); ++nodeIter) {

	}

	// write to file

}

void nodelist::writePOV() {
	std::stringstream fileNameStream;
	fileNameStream << "./data/rods_" << name << ".pov";
	std::string outfileName = fileNameStream.str();

	std::vector<std::string> data;

	for (std::vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		std::stringstream stream;

		// if point is without bounary...
		if ((*nodeIter)->getPosition().min() < -5
				|| (*nodeIter)->getPosition().max() > 5) {
			continue;
		}

		// sphere at joint
		stream << "sphere{"
				<< (*nodeIter)->getPosition().toString("<", ",", ">") << ",r}"
				<< std::endl;

		for (std::vector<node*>::iterator neighIter =
				(*nodeIter)->getNeighbours()->begin();
				neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

			// get target link
			coordinate linkTarget(std::numeric_limits<double>::infinity(),
					std::numeric_limits<double>::infinity(),
					std::numeric_limits<double>::infinity());
			if (periodic) {

				double distSqr = std::numeric_limits<double>::infinity();
				coordinate diff;

				std::vector<coordinate> shifters = getShifters();
				for (std::vector<coordinate>::iterator shifter =
						shifters.begin(); shifter != shifters.end();
						++shifter) {
					diff = (((*neighIter)->getPosition()) + (*shifter))
							- ((*nodeIter)->getPosition());
					if (diff.lengthSqr() < distSqr) {
						distSqr = diff.lengthSqr();
						linkTarget =
								(((*neighIter)->getPosition()) + (*shifter));
					}
				}
			} else {
				linkTarget = (*neighIter)->getPosition();
			}

			// sphere at joint over the boundary
			stream << "sphere{" << linkTarget.toString("<", ",", ">") << ",r}"
					<< std::endl;

			// cylinder from a to b TODO: den zurück nicht...
			stream << "cylinder{"
					<< (*nodeIter)->getPosition().toString("<", ",", ">") << ","
					<< linkTarget.toString("<", ",", ">") << ",r}" << std::endl;

		}
		data.push_back(stream.str());
	}

	// write to outfile
	std::ofstream outfile;
	outfile.open(outfileName.c_str());
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	std::cout << "Pattern written to file '" << outfileName
			<< "', ready to be rendered by POV-Ray." << std::endl;

}

void nodelist::writeMEEP() {
	std::stringstream fileNameStream;
	fileNameStream << "./data/meep-dielectric_" << name << ".ctl";
	std::string outfileName = fileNameStream.str();

	std::vector<std::string> data;

	// open list
	data.push_back("(list\n");

	// Iteration over nodes
	for (std::vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		std::stringstream stream;
		stream.precision(4);

		// Cylinders for rods
		//double radius = 0.2685;
		stream << "(make sphere (material polymer) (center "
				<< (*nodeIter)->getPosition().toString("", " ", "")
				<< ") (radius rad)) ";

		for (std::vector<node*>::iterator neighIter =
				(*nodeIter)->getNeighbours()->begin();
				neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

			// TODO: as extra function, as it is used quite often...
			// get target link
			coordinate linkTarget(std::numeric_limits<double>::infinity(),
					std::numeric_limits<double>::infinity(),
					std::numeric_limits<double>::infinity());
			if (periodic) {

				double distSqr = std::numeric_limits<double>::infinity();
				coordinate diff;

				std::vector<coordinate> shifters = getShifters();
				for (std::vector<coordinate>::iterator shifter =
						shifters.begin(); shifter != shifters.end();
						++shifter) {
					diff = (((*neighIter)->getPosition()) + (*shifter))
							- ((*nodeIter)->getPosition());
					if (diff.lengthSqr() < distSqr) {
						distSqr = diff.lengthSqr();
						linkTarget =
								(((*neighIter)->getPosition()) + (*shifter));
					}
				}

			} else {
				linkTarget = (*neighIter)->getPosition();
			}
			// e2 (width), orthogonal to e3 (height) and e1 (length = connecting vector) => cross product
			coordinate axis = (*nodeIter)->getPosition() - linkTarget;

			// length = length of connecvting vector
			double height = axis.length();

			// center: middpoint of points
			std::string center =
					(((*nodeIter)->getPosition() + linkTarget) / 2).toString("",
							" ", "");

			// block from a to b TODO: den zurück nicht...
			stream << "(make cylinder (material polymer) (center " << center
					<< ") (radius rad) (height " << height << ") (axis "
					<< axis.toString("", " ", "") << ")) ";
		}

		// Ellipsoid Approximation (apparently too big ctl file for meep...)
		/*

		 // assuiming a writing spheroid with dimensions w x w x h
		 const double w = 0.2;
		 const double h = 0.6;
		 // ellipsoid at joint
		 stream << "(make ellipsoid (material polymer) (center "
		 << (*nodeIter)->getPosition().toString("", " ", "")
		 << ") (size " << w << " " << w << " " << h << "))";// << endl;
		 // number of blocks for ellisoid approximation
		 int numBlocks = 5;
		 for (vector<node*>::iterator neighIter =
		 (*nodeIter)->getNeighbours()->begin();
		 neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

		 // cylinder only if not connected over boundary
		 if ((!periodic) || ((*nodeIter)->euklidian((*neighIter))
		 < getMaxFeatureSize())) { {

		 // e2 (width), orthogonal to e3 (height) and e1 (length = connecting vector) => cross product
		 coordinate e1 = (*nodeIter)->getPosition()
		 - (*neighIter)->getPosition();
		 coordinate e3(0, 0, 1);
		 coordinate e2 = coordinate::cpr(e1, e3);

		 // length = length of connecvting vector
		 double l = e1.length();

		 // center: middpoint of points
		 string center = (((*nodeIter)->getPosition()
		 + (*neighIter)->getPosition()) / 2).toString("", " ",
		 "");

		 // ellipsoid approximation
		 for (double theta = 0; theta <= M_PI / 2;
		 theta += M_PI / (2 * numBlocks)) {
		 // dimensions of blocks used to approximnate ellipse
		 double wApp = 2 * r(theta, w, h) * cos(theta);
		 double hApp = 2 * r(theta, w, h) * sin(theta);

		 // block from a to b TODO: den zurück nicht...
		 stream << "(make block (material polymer) (center "
		 << center << ") (size " << l << " " << wApp << " "
		 << hApp << ") (e1 " << e1.toString("", " ", "")
		 << ") (e2 " << e2.toString("", " ", "") << ") (e3 "
		 << e3.toString("", " ", "") << "))";// << endl;
		 }

		 }
		 }*/
		stream << std::endl;

		data.push_back(stream.str());
	}

// close list
	data.push_back(")\n");

// write to outfile
	std::ofstream outfile;
	outfile.open(outfileName.c_str());
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	std::cout << "Pattern written as " << outfileName << ", ready to be MEEPed."
			<< std::endl;
}

void nodelist::writeMPB() {

	const char outfileName[] = "./data/mpb-dielectric.ctl";

	std::vector<std::string> data;

	// open list
	data.push_back("(list\n");

	// Iteration over nodes
	for (std::vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		std::stringstream stream;
		stream.precision(4);

		// Cylinders for rods
		//double radius = 0.2685;
		stream
				<< (*nodeIter)->getPosition().toString(
						"(make sphere (material polymer) (center (c->l ", " ",
						")) (radius rad)) ");

		for (std::vector<node*>::iterator neighIter =
				(*nodeIter)->getNeighbours()->begin();
				neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

			// cylinder only if not connected over boundary
			// !a || (a && b) == !a || b (a - periodic, b - link goes over edge)
			if ((!periodic)
					|| ((*nodeIter)->euklidian((*neighIter))
							< getMaxFeatureSize())) {

				// e2 (width), orthogonal to e3 (height) and e1 (length = connecting vector) => cross product
				coordinate axis = (*nodeIter)->getPosition()
						- (*neighIter)->getPosition();

				// length = length of connecvting vector
				double height = axis.length();

				// center: middpoint of points
				std::string center = (((*nodeIter)->getPosition()
						+ (*neighIter)->getPosition()) / 2).toString("", " ",
						"");

				// cylinder from a to b TODO: den zurück nicht...
				stream << "(make cylinder (material polymer) (center (c->l "
						<< center << ")) (radius rad) (height " << height
						<< ") (axis " << axis.toString("(c->l ", " ", ")")
						<< ")) ";

			}
		}

		// Ellipsoid Approximation (apparently too big ctl file for meep...)
		/*

		 // assuiming a writing spheroid with dimensions w x w x h
		 const double w = 0.2;
		 const double h = 0.6;
		 // ellipsoid at joint
		 stream << "(make ellipsoid (material polymer) (center "
		 << (*nodeIter)->getPosition().toString("", " ", "")
		 << ") (size " << w << " " << w << " " << h << "))";// << endl;
		 // number of blocks for ellisoid approximation
		 int numBlocks = 5;
		 for (vector<node*>::iterator neighIter =
		 (*nodeIter)->getNeighbours()->begin();
		 neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

		 // cylinder only if not connected over boundary
		 if ((!periodic) || ((*nodeIter)->euklidian((*neighIter))
		 < getMaxFeatureSize())) { {

		 // e2 (width), orthogonal to e3 (height) and e1 (length = connecting vector) => cross product
		 coordinate e1 = (*nodeIter)->getPosition()
		 - (*neighIter)->getPosition();
		 coordinate e3(0, 0, 1);
		 coordinate e2 = coordinate::cpr(e1, e3);

		 // length = length of connecvting vector
		 double l = e1.length();

		 // center: middpoint of points
		 string center = (((*nodeIter)->getPosition()
		 + (*neighIter)->getPosition()) / 2).toString("", " ",
		 "");

		 // ellipsoid approximation
		 for (double theta = 0; theta <= M_PI / 2;
		 theta += M_PI / (2 * numBlocks)) {
		 // dimensions of blocks used to approximnate ellipse
		 double wApp = 2 * r(theta, w, h) * cos(theta);
		 double hApp = 2 * r(theta, w, h) * sin(theta);

		 // block from a to b TODO: den zurück nicht...
		 stream << "(make block (material polymer) (center "
		 << center << ") (size " << l << " " << wApp << " "
		 << hApp << ") (e1 " << e1.toString("", " ", "")
		 << ") (e2 " << e2.toString("", " ", "") << ") (e3 "
		 << e3.toString("", " ", "") << "))";// << endl;
		 }

		 }
		 }*/
		stream << std::endl;

		data.push_back(stream.str());
	}

// close list
	data.push_back(")\n");

// write to outfile
	std::ofstream outfile;
	outfile.open(outfileName);
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	std::cout << "Pattern written as " << outfileName << ", ready to be MPBed."
			<< std::endl;
}

std::vector<node*>::iterator nodelist::begin() {
	return list.begin();
}

std::vector<node*>::iterator nodelist::end() {
	return list.end();
}

void nodelist::setEdgenodes(double distance) {
	for (std::vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		(*n)->setEdgenode(distance);
	}
}

node& nodelist::operator [](int i) {
	return *list[i];
}
/*
 const node& nodelist::operator [](const int i) const {
 return &list[i];
 }
 */
int nodelist::size() {
	return list.size();
}

std::vector<std::vector<double> > nodelist::getGnuplotMatrix() {
	std::vector<double> emptyLine;
	for (unsigned int i = 0; i < 3; i++) {
		emptyLine.push_back(nan(""));
	}

	std::vector<std::vector<double> > data;
	// nodesiteration
	for (std::vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		for (unsigned int n = 0; n < (*nodeIter)->countNeighbours(); n++) {

			coordinate shifter = (*nodeIter)->getNeighbourShifters()->at(n);
			// do not plot links over the boundary
			if ((*nodeIter)->getNeighbourShifters()->at(n)
					== coordinate(0, 0, 0)) {
				data.push_back(*(*nodeIter)->getPosition().getVector());
				data.push_back(
						*((*nodeIter)->getNeighbours()->at(n)->getPosition()
								- shifter).getVector());

				data.push_back(emptyLine);
			}
		}

		// neighboursiteration
		/*		for (std::vector<node*>::iterator neighIter =
		 (*nodeIter)->getNeighbours()->begin();
		 neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

		 // get target link
		 coordinate linkTarget(std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity());
		 if (periodic) {

		 double distSqr = std::numeric_limits<double>::infinity();
		 coordinate diff;

		 std::vector<coordinate> shifters = getShifters();
		 for (std::vector<coordinate>::iterator shifter =
		 shifters.begin(); shifter != shifters.end();
		 ++shifter) {
		 diff = (((*neighIter)->getPosition()) + (*shifter))
		 - ((*nodeIter)->getPosition());
		 if (diff.lengthSqr() < distSqr) {
		 distSqr = diff.lengthSqr();
		 linkTarget =
		 (((*neighIter)->getPosition()) + (*shifter));
		 }
		 }

		 } else {
		 linkTarget = (*neighIter)->getPosition();
		 }

		 data.push_back(*(*nodeIter)->getPosition().getVector());
		 data.push_back(*linkTarget.getVector());

		 data.push_back(emptyLine);

		 }*/
	}
	return data;
}

std::vector<std::vector<double> > nodelist::getEdgelinksGnuplotMatrix() {
	std::vector<double> emptyLine;
	for (unsigned int i = 0; i < 3; i++) {
		emptyLine.push_back(nan(""));
	}

	std::vector<std::vector<double> > data;
	// nodesiteration
	for (std::vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		// neighboursiteration TODO: still an error somewhere... possibly in the shifter assigning
		for (unsigned int n = 0; n < (*nodeIter)->countNeighbours(); n++) {

			coordinate shifter = (*nodeIter)->getNeighbourShifters()->at(n);
			// its only a link over the boundary if there's not the 0-shifter...
			if ((*nodeIter)->getNeighbourShifters()->at(n)
					!= coordinate(0, 0, 0)) {
				data.push_back(*(*nodeIter)->getPosition().getVector());
				data.push_back(
						*((*nodeIter)->getNeighbours()->at(n)->getPosition()
								- shifter).getVector());

				data.push_back(emptyLine);
			}
			/*

			 for (std::vector<node*>::iterator neighIter =
			 (*nodeIter)->getNeighbours()->begin();
			 neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

			 // get target link
			 coordinate linkTarget(std::numeric_limits<double>::infinity(),
			 std::numeric_limits<double>::infinity(),
			 std::numeric_limits<double>::infinity());

			 double distSqr = std::numeric_limits<double>::infinity();
			 coordinate diff;

			 coordinate usedShifter = coordinate(0, 0, 0);

			 // shift the neighbour with the periodic boundary conditions until its closest to the node
			 std::vector<coordinate> shifters = getShifters();
			 for (std::vector<coordinate>::iterator shifter = shifters.begin();
			 shifter != shifters.end(); ++shifter) {
			 diff = (((*neighIter)->getPosition()) + (*shifter))
			 - ((*nodeIter)->getPosition());
			 if (diff.lengthSqr() < distSqr) {
			 distSqr = diff.lengthSqr();
			 linkTarget = (((*neighIter)->getPosition()) + (*shifter));
			 usedShifter = (*shifter);
			 }
			 }

			 // exclude the zero shifter
			 if (usedShifter != coordinate(0, 0, 0)) {
			 data.push_back(*(*nodeIter)->getPosition().getVector());
			 data.push_back(*linkTarget.getVector());

			 data.push_back(emptyLine);
			 }
			 */
		}
	}
	return data;
}

void nodelist::deleteEntries() {
	for (std::vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		delete (*n);
	}
}

