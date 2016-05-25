/*
 * nodelist.cpp
 *
 *  Created on: 24.07.2015
 *      Author: Lukas Siedentop
 */

#include "nodelist.hpp"

using namespace std;

nodelist::nodelist() :
		periodic(0), name(""), min(coordinate()), max(coordinate()) {
}

nodelist::nodelist(bool periodicity, string n) :
		periodic(periodicity), name(n), min(coordinate()), max(coordinate()) {
}

nodelist::nodelist(int pattern, bool periodicity) :
		periodic(periodicity), min(coordinate(-5, -5, -5)), max(
				coordinate(5, 5, 5)) {

	switch (pattern) {
	case 1: { // poisson pattern
		cout << "Generate poisson pattern with 1000 points in 10^3 cubicle..."
				<< endl;

		name = "poissonpattern";

		// seed for random numbers generator
		srand(time(NULL));
		for (unsigned n = 0; n < 1000; n++) {
			node* tmp = new node(this,
					(coordinate(3) * 10) - coordinate(5, 5, 5));
			tmp->setEdgenode(1);
			list.push_back(tmp);
		}
		cout << "Poisson pattern";
		break;
	}
	case 2: { // diamond point pattern ((c) Dirk)
		cout
				<< "Generate diamond point pattern with 1000 points in 10^3 cubicle..."
				<< endl;

		name = "diamondpattern";

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
		for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
			(*n)->setEdgenode(1);
		}

		cout << "Diamond point pattern";
		break;
	}
	}
	setNeighbours();

	cout << " generated. Statistics:" << endl << listStats() << endl;
}

void nodelist::setNeighbours() {
	// get four next points
	for (unsigned int neighs = 0; neighs < 4; neighs++) {
		// nodesiteration
		for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
			node* neigh = NULL;
			if ((*n)->countNeighbours() < 4) {

				// find closest neighbour with less than four neighbours
				double distSqr = std::numeric_limits<double>::infinity();
				// 2nd nodeiteration
				for (vector<node*>::iterator nextn = list.begin();
						nextn != list.end(); ++nextn) {

					if ((*n)->isNeighbour((*nextn)) || n == nextn
							|| (*nextn)->countNeighbours() >= 4) {
						continue;
					}

					coordinate diff;

					// compare the difference vectors with an added shifter for the periodic case
					if (periodic) {
						vector<coordinate> shifters = getShifters();
						for (vector<coordinate>::iterator shifter =
								shifters.begin(); shifter != shifters.end();
								++shifter) {
							diff = (((*nextn)->getPosition()) + (*shifter))
									- ((*n)->getPosition());
							if (diff.lengthSqr() < distSqr) {
								distSqr = diff.lengthSqr();
								neigh = (*nextn);
							}
						}
					} else {
						diff = ((*nextn)->getPosition())
								- ((*n)->getPosition());
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
		}
	}
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
	for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
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

int nodelist::pointsInside(const vector<coordinate>& points,
		const coordinate& mid, const double r) const {
	int ctr = 0;

	// nodesiterations
	for (vector<coordinate>::const_iterator point = points.begin();
			point != points.end(); ++point) {
		// bounding box of the sphere
		if ((fabs((*point)[0] - mid[0]) < r) && (fabs((*point)[1] - mid[1]) < r)
				&& (fabs((*point)[2] - mid[2]) < r)) {
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
	for (vector<node*>::iterator it = list.begin(); it != list.end(); ++it) {
		(*it)->shift(shifter);
	}
}

void nodelist::scaleList(double a) {
	// correct bounding box
	min *= a;
	max *= a;

	// scale entries
	for (vector<node*>::iterator it = list.begin(); it != list.end(); ++it) {
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
	for (vector<node*>::iterator it = list.begin(); it != list.end(); ++it) {
		(*it)->scaleAnisotropic(ax, ay, az);
	}
}

/* unneccesary
 bool nodelist::isPeriodic() const {
 return periodic;
 }
 */

string nodelist::getName() {
	return name;
}

node* nodelist::add(double x, double y, double z) {

	// check if point at given position exist
	for (vector<node*>::iterator it = list.begin(); it != list.end(); ++it) {
		if (periodic) {
			coordinate pos = ((*it)->getPosition());
			if ((std::min(fabs(x - pos[0]),
					(fabs(fabs(x - pos[0]) - getLengths()[0]))) < 0.000001)
					&& (std::min(fabs(y - pos[1]),
							(fabs(fabs(y - pos[1]) - getLengths()[1])))
							< 0.000001)
					&& (std::min(fabs(z - pos[2]),
							(fabs(fabs(z - pos[2]) - getLengths()[2])))
							< 0.000001)) {
				return *it;
			}
		}
		if (!periodic && ((*it)->getPosition() == coordinate(x, y, z))) {
			return *it;
		}
	}

	coordinate newPos = coordinate(x, y, z);

	// ensure that point lies in between boundaries if periodic
	if (periodic) {
		vector<coordinate> shifters = getShifters();
		for (vector<coordinate>::iterator shifter = shifters.begin();
				shifter != shifters.end(); ++shifter) {
			// as bounding box
			if ((fabs((*shifter)[0] + newPos[0]) < getLengths()[0] / 2)
					&& (fabs((*shifter)[1] + newPos[1]) < getLengths()[1] / 2)
					&& (fabs((*shifter)[2] + newPos[2]) < getLengths()[2] / 2)) {
				newPos += (*shifter);
			}
		}
	}

	// if not, build and add new node
	node* n = new node(this, newPos);
	list.push_back(n);
	return n;
}
/* unnecessary
 coordinate nodelist::getMins() {
 return min;
 }

 coordinate nodelist::getMaxs() {
 return max;
 }*/

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
	cout << "Number of elements: " << list.size() << endl;

	int i = 1;

	// nodesiteration
	for (vector<node*>::iterator it = list.begin(); it != list.end(); ++it) {
		cout << char(9) << i << "-th element: " << ((*it)->getPosition())
				<< " with " << (*it)->countNeighbours() << " neighbours at: ";
		// neighbouriteration
		for (vector<class node *>::iterator neighIt =
				(*it)->getNeighbours()->begin();
				neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
			cout << ((*neighIt)->getPosition()) << " ";
		}
		cout << endl;

		i++;
	}
}

vector<coordinate> nodelist::getShifters() {
	vector<coordinate> shifters;

	// boxlength
	double lx = getLengths()[0], ly = getLengths()[1], lz = getLengths()[2];

	// go through all combinations
	for (int ix = -1; ix < 2; ix++) {
		for (int iy = -1; iy < 2; iy++) {
			for (int iz = -1; iz < 2; iz++) {
				// add the shifting vector
				shifters.push_back(coordinate(ix * lx, iy * ly, iz * lz));
			}
		}
	}
	return shifters;
}

vector<coordinate> nodelist::getShifted(coordinate mid, double halfExtend) {
	vector<coordinate> shifters;

	// boxlength
	double lx = getLengths()[0], ly = getLengths()[1], lz = getLengths()[2]; // 5

	// go through all combinations, but only if one edge of the given box is not within the pattern...
	for (int ix = -1; ix < 2; ix++) {
		if ((fabs(mid[0] + (ix * halfExtend)) >= fabs(ix * (lx / 2.0)))) {
			for (int iy = -1; iy < 2; iy++) {
				if ((fabs(mid[1] + (iy * halfExtend)) >= fabs(iy * (ly / 2.0)))) {
					for (int iz = -1; iz < 2; iz++) {
						if ((fabs(mid[2] + (iz * halfExtend))
								>= fabs(iz * (lz / 2.0)))) {
							// ... add the shifting vector (Shift sphere contrary to pattern!)
							shifters.push_back(
									mid
											- coordinate(ix * lx, iy * ly,
													iz * lz));
						}
					}
				}
			}
		}
	}
	return shifters;
}

string nodelist::listStats(const string commentDelimeter) {
	stringstream stream;
	stream << commentDelimeter << (periodic ? "Periodic" : "Non-periodic")
			<< endl;
	if (!periodic) {
		stream << commentDelimeter << "Number of edgenodes: "
				<< countEdgenodes() << endl;
	}


	stream << commentDelimeter << "Number of nodes: " << list.size() << endl;
	stream << commentDelimeter << "Bounding box: " << min << ", " << max
			<< endl;
	stream << commentDelimeter << "Characteristic length: " << stats(lengthDistribution())[1] << endl;
	stream << commentDelimeter << "Mid of box: " << getMid() << endl;
	stream << commentDelimeter << "Volume: " << getVolume() << endl;
	stream << commentDelimeter << "Density of points: " << getDensity() << endl;

	return stream.str();
}

double nodelist::normalize() {
	// stats(lengthDistribution())[1]

	// only for agapornis! TODO: function to set scaling for each direction...
	double resX = 9.067;// nm/px
	double resY = 9.067;// nm/px
	double resZ = 25;//nm/px

	scaleListAnisotropic(resX, resY, resZ);

	// set density to 1
	setDensity(1);

	// midpoint
	shiftList(getMid() * -1);

	double factor = pow(getDensity(), 1. / 3.);
	return factor;
}

vector<double> nodelist::neighbourDistribution() {
	vector<double> data;

	int counter;

	// nodesiteration
	for (vector<node*>::iterator it = list.begin(); it != list.end(); ++it) {
		if (!(*it)->isEdgenode()) {
			counter = 0;
			// neighboursiteration
			for (vector<node*>::iterator neighIt =
					(*it)->getNeighbours()->begin();
					neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
				counter++;
			}

			data.push_back(counter);
		}
	}

	return data;
}

vector<double> nodelist::lengthDistribution() {
	vector<double> data;

	// nodesiteration
	for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		// exclude distances between edgendoes
		if (!(*n)->isEdgenode()) {
			// neighboursiteration
			for (vector<node*>::iterator nn = (*n)->getNeighbours()->begin();
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
	vector<double> halfData;
	for (unsigned int i = 0; i < data.size(); i += 2) {
		halfData.push_back(data[i]);
	}

	return halfData;
}

vector<double> nodelist::angleDistribution() {
	vector<double> data;

	// nodesiteration
	for (vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		if (!(*nodeIter)->isEdgenode()) {
			// neighboursiteration
			for (vector<node*>::iterator neighIter1 =
					(*nodeIter)->getNeighbours()->begin();
					neighIter1 != (*nodeIter)->getNeighbours()->end();
					++neighIter1) {

				// neighbours of neighboursiteration
				for (vector<node*>::iterator neighIter2 =
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
	vector<double> halfData;
	for (unsigned int i = 0; i < data.size(); i += 2) {
		halfData.push_back(data[i]);
	}

	return halfData;
}

vector<vector<double> > nodelist::hyperuniformity(unsigned int nr,
		unsigned int n) {
	// radiusincrement + maximal radius
	double rMax = periodic ? getLengths().min() : getLengths().min() / 3, dr =
			rMax / nr;

	// datastructure to save the number of points in each sphere
	vector<vector<double> > data;
	data.resize(nr);
	for (unsigned int i = 0; i < nr; ++i) {
		data[i].resize(n);
	}
	clock_t t;
	t = clock();

	// seed for generation of random numbers
	srand(time(NULL));

	// midpoint of sphere
	coordinate mid;

	if (periodic) {
		// precalculate extended pattern
		vector<coordinate> shifters = getShifters();
		vector<coordinate> extendedPattern;
		for (vector<coordinate>::iterator shifter = shifters.begin();
				shifter != shifters.end(); ++shifter) {
			for (vector<node*>::iterator n = list.begin(); n != list.end();
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
			for (unsigned int j = 0; j < nr; j++) {
				data[j][i] = pointsInside(extendedPattern, mid, j * dr);
			}

			++show_progress;
		}
	} else {
		// get patterns coordinates
		vector<coordinate> pattern;
		for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
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
	cout << "Took " << t << " Clicks (" << (((float) t) / CLOCKS_PER_SEC)
			<< "s)" << endl;

	//cout << "Radius" << char(9) << "Volumen" << char(9) << "Erwartungswert" << char(9) << "Verhältnis" << endl;
	// Calculate the expected value (of numbers of points within sphere) of each radius. Should roughly be equal to the numberdensity times volume of the sphere.
	vector<double> expectedValue;
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
	vector<vector<double> > variance;
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
	stringstream fileNameStream;
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
	for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		stringstream stream;
		stream << ((*n)->getPosition()).toString("", "\t", "", 10) << endl;
		data.push_back(stream.str());
	}

	// write to outfile
	ofstream outfile;
	outfile.open(outfileName.c_str());
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	cout << "Pattern written to file '" << outfileName << "'. Enjoy!" << endl;
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
	vector<coordinate> loop;
	for (vector<node*>::iterator nodeIter = extendedList.begin();
			nodeIter != extendedList.end(); ++nodeIter) {

	}

	// write to file

}

void nodelist::writePOV() {
	stringstream fileNameStream;
	fileNameStream << "./data/rods_" << name << ".pov";
	string outfileName = fileNameStream.str();

	vector<string> data;

	for (vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		stringstream stream;

		// if point is without bounary...
		if ((*nodeIter)->getPosition().min() < -5 || (*nodeIter)->getPosition().max() > 5) {
			continue;
		}

		// sphere at joint
		stream << "sphere{"
				<< (*nodeIter)->getPosition().toString("<", ",", ">") << ",r}"
				<< endl;

		for (vector<node*>::iterator neighIter =
				(*nodeIter)->getNeighbours()->begin();
				neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

			// get target link
			coordinate linkTarget(std::numeric_limits<double>::infinity(),
					std::numeric_limits<double>::infinity(),
					std::numeric_limits<double>::infinity());
			if (periodic) {

				double distSqr = std::numeric_limits<double>::infinity();
				coordinate diff;

				vector<coordinate> shifters = getShifters();
				for (vector<coordinate>::iterator shifter = shifters.begin();
						shifter != shifters.end(); ++shifter) {
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
			stream << "sphere{"
					<< linkTarget.toString("<", ",", ">") << ",r}"
					<< endl;

			// cylinder from a to b TODO: den zurück nicht...
			stream << "cylinder{"
					<< (*nodeIter)->getPosition().toString("<", ",", ">") << ","
					<< linkTarget.toString("<", ",", ">") << ",r}" << endl;

		}
		data.push_back(stream.str());
	}

	// write to outfile
	ofstream outfile;
	outfile.open(outfileName.c_str());
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	cout << "Pattern written to file '" << outfileName
			<< "', ready to be rendered by POV-Ray." << endl;

}

void nodelist::writeMEEP() {
	stringstream fileNameStream;
	fileNameStream << "./data/meep-dielectric_" << name << ".ctl";
	string outfileName = fileNameStream.str();

	vector<string> data;

	// open list
	data.push_back("(list\n");

	// Iteration over nodes
	for (vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		stringstream stream;
		stream.precision(4);

		// Cylinders for rods
		//double radius = 0.2685;
		stream << "(make sphere (material polymer) (center "
				<< (*nodeIter)->getPosition().toString("", " ", "")
				<< ") (radius rad)) ";

		for (vector<node*>::iterator neighIter =
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

				vector<coordinate> shifters = getShifters();
				for (vector<coordinate>::iterator shifter = shifters.begin();
						shifter != shifters.end(); ++shifter) {
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
			string center =
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
		stream << endl;

		data.push_back(stream.str());
	}

// close list
	data.push_back(")\n");

// write to outfile
	ofstream outfile;
	outfile.open(outfileName.c_str());
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	cout << "Pattern written as " << outfileName << ", ready to be MEEPed."
			<< endl;
}

void nodelist::writeMPB() {

	const char outfileName[] = "./data/mpb-dielectric.ctl";

	vector<string> data;

	// open list
	data.push_back("(list\n");

	// Iteration over nodes
	for (vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		stringstream stream;
		stream.precision(4);

		// Cylinders for rods
		//double radius = 0.2685;
		stream
				<< (*nodeIter)->getPosition().toString(
						"(make sphere (material polymer) (center (c->l ", " ",
						")) (radius rad)) ");

		for (vector<node*>::iterator neighIter =
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
				string center = (((*nodeIter)->getPosition()
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
		stream << endl;

		data.push_back(stream.str());
	}

// close list
	data.push_back(")\n");

// write to outfile
	ofstream outfile;
	outfile.open(outfileName);
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	cout << "Pattern written as " << outfileName << ", ready to be MPBed."
			<< endl;
}

std::vector<node*>::iterator nodelist::begin() {
	return list.begin();
}

std::vector<node*>::iterator nodelist::end() {
	return list.end();
}

void nodelist::setEdgenodes(double distance) {
	for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
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

vector<vector<double> > nodelist::getGnuplotMatrix() {
	vector<double> emptyLine;
	for (unsigned int i = 0; i < 3; i++) {
		emptyLine.push_back(nan(""));
	}

	vector<vector<double> > data;
	// nodesiteration
	for (vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		// neighboursiteration
		for (vector<node*>::iterator neighIter =
				(*nodeIter)->getNeighbours()->begin();
				neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

			// get target link
			coordinate linkTarget(std::numeric_limits<double>::infinity(),
					std::numeric_limits<double>::infinity(),
					std::numeric_limits<double>::infinity());
			if (periodic) {

				double distSqr = std::numeric_limits<double>::infinity();
				coordinate diff;

				vector<coordinate> shifters = getShifters();
				for (vector<coordinate>::iterator shifter = shifters.begin();
						shifter != shifters.end(); ++shifter) {
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

		}
	}
	return data;
}

void nodelist::deleteEntries() {
	for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		delete (*n);
	}
}

