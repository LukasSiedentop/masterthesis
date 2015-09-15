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
	case 1: { // random point pattern
		cout
				<< "Generate random point pattern with 1000 points in 10^3 cubicle..."
				<< endl;

		name = "random point pattern";

		// seed for random numbers generator
		srand(time(NULL));
		for (unsigned n = 0; n < 1000; n++) {
			node* tmp = new node(this,
					(coordinate(3) * 10) - coordinate(5, 5, 5));
			tmp->setEdgenode(1);
			list.push_back(tmp);
		}
		cout << "Random point pattern";
		break;
	}
	case 2: { // diamond point pattern ((c) Dirk)
		cout
				<< "Generate diamond point pattern with 1000 points in 10^3 cubicle..."
				<< endl;

		name = "diamond point pattern";

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
							diff = ((*nextn)->getPosition())
									- ((*n)->getPosition()) + (*shifter);
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

bool nodelist::isPeriodic() const {
	return periodic;
}

void nodelist::setDensity(double density) {
	scaleList(pow((getDensity() / density), 1. / 3.));
}

string nodelist::getName() {
	return name;
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

	// if not build and add new node
	node* n = new node(this, coordinate(x, y, z));
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

int nodelist::pointsInside(coordinate mid, double r, double rSqr) {
	int ctr = 0;

	// nodesiterations
	for (vector<node*>::iterator it = list.begin(); it != list.end(); ++it) {
		// bounding box of the sphere
		if ((fabs(((*it)->getPosition())[0] - mid[0]) < r)
				&& (fabs(((*it)->getPosition())[1] - mid[1]) < r)
				&& (fabs(((*it)->getPosition())[2] - mid[2]) < r)) {
			// sphere itself
			if ((((*it)->getPosition()) - mid).lengthSqr() < rSqr) {
				ctr++;
			}
		}
	}
	return ctr;
}

int nodelist::pointsInsidePeriodic(coordinate mid, double r) {
	int ctr = 0;

	vector<coordinate> shifted = getShifted(mid, r);

	// distance comparison of sphere shifted in every dimension
	for (unsigned int s = 0; s < shifted.size(); s++) {
		// nodesiteration
		for (vector<node*>::iterator it = list.begin(); it != list.end();
				++it) {
			// bounding box of sphere
			if ((fabs(((*it)->getPosition())[0] + shifted[s][0]) < r)
					&& (fabs(((*it)->getPosition())[1] + shifted[s][1]) < r)
					&& (fabs(((*it)->getPosition())[2] + shifted[s][2]) < r)) {

				// Cubicles are not stricly convex, unfortunately... (Local density fluctuations, hyperuniformity, and order metrics. Salvatore Torquato et. al.)
				//ctr++;

				// is the point within the sphere?
				if (((((*it)->getPosition()) + shifted[s]).lengthSqr())
						< (r * r)) {
					ctr++;
				}
			}
		}
	}
	return ctr;
}

string nodelist::listStats(const char commentDelimeter[]) {
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
	stream << commentDelimeter << "Mid of box: " << getMid() << endl;
	stream << commentDelimeter << "Volume: " << getVolume() << endl;
	stream << commentDelimeter << "Density of points: " << getDensity() << endl;

	return stream.str();
}

double nodelist::normalize() {
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
					// dont calculate angles with itself
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

vector<vector<double> > nodelist::hyperuniformity(unsigned int nr, unsigned int n) {
	// number of radii
	//unsigned int nr = 50;
	// radiusincrement + maximal radius
	double rMax = periodic ? getLengths().min() : getLengths().min() / 3, dr =
			rMax / nr;

	// number of spheres/iterations
	//unsigned int n = 10000;

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
	coordinate mid;
	if (periodic) {
		boost::progress_display show_progress(n);

		// iteration over n spheres
		for (unsigned int i = 0; i < n; i++) {
			// Choose the center of the sphere such that it is somewhere within the whole pattern.
			mid = (coordinate(3) - 0.5) * getLengths(); // faster, TODO not tested
					/*mid = coordinate(3); old
					 mid *= getLengths();
					 mid -= getLengths() / 2;*/
			// Without this, the Pattern is expected to be centered around (0,0,0).
			//mid += getMid();

			// iteration over radius: count points within sphere
			for (unsigned int j = 0; j < nr; j++) {
				data[j][i] = pointsInsidePeriodic(mid, j * dr);
			}

			++show_progress;
		}
	} else {
		boost::progress_display show_progress(nr);

		double r, rSqr;
		// iteration over radius
		for (unsigned int j = 0; j < nr; j++) {
			r = j * dr;
			rSqr = r * r;
			// iteration over n spheres
			for (unsigned int i = 0; i < n; i++) {
				// Choose the center of the sphere such that it is somewhere within the whole pattern.
				mid = (coordinate(3) - 0.5) * (getLengths() - 2 * r);

				data[j][i] = pointsInside(mid, r, rSqr);
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

void nodelist::writePOV() {

	const char outfileName[] = "./data/staebe.pov";

	vector<string> data;

	for (vector<node*>::iterator nodeIter = list.begin();
			nodeIter != list.end(); ++nodeIter) {
		stringstream stream;

		// sphere at joint
		stream << "sphere{"
				<< (*nodeIter)->getPosition().toString("<", ",", ">") << ",r}"
				<< endl;

		for (vector<node*>::iterator neighIter =
				(*nodeIter)->getNeighbours()->begin();
				neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

			// cylinder only if not connected over boundary
			if (isPeriodic()
					&& (*nodeIter)->euklidian((*neighIter))
							< getMaxFeatureSize()) {
				// cylinder from a to b TODO: den zurück nicht...
				stream << "cylinder{"
						<< (*nodeIter)->getPosition().toString("<", ",", ">")
						<< ","
						<< (*neighIter)->getPosition().toString("<", ",", ">")
						<< ",r}" << endl;

			} else if (!isPeriodic()) {
				// cylinder from a to b TODO: den zurück nicht...
				stream << "cylinder{"
						<< (*nodeIter)->getPosition().toString("<", ",", ">")
						<< ","
						<< (*neighIter)->getPosition().toString("<", ",", ">")
						<< ",r}" << endl;
			}
		}
		data.push_back(stream.str());
	}

	// write to outfile
	ofstream outfile;
	outfile.open(outfileName);
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}
	outfile.close();

	cout << "Pattern written in " << outfileName
			<< ", ready to be rendered by POV-Ray." << endl;

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

			//  print link only if both are not connected over boundary
			if (periodic
					&& (*nodeIter)->euklidian(*neighIter)
							< getMaxFeatureSize()) {

				data.push_back(*(*nodeIter)->getPosition().getVector());
				data.push_back(*(*neighIter)->getPosition().getVector());

				data.push_back(emptyLine);

			} else if (!periodic) {
				//if (!(*neighIter)->isEdgenode()) {
				//if (!(*nodeIter)->isEdgenode()) {
				data.push_back(*(*nodeIter)->getPosition().getVector());
				data.push_back(*(*neighIter)->getPosition().getVector());

				data.push_back(emptyLine);

				//}
				//}
			}
		}
	}
	return data;
}

void nodelist::deleteEntries() {
	for (vector<node*>::iterator n = list.begin(); n != list.end(); ++n) {
		delete (*n);
	}
}

