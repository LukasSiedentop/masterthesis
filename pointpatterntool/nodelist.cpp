/*
 * nodelist.cpp
 *
 *  Created on: 24.07.2015
 *      Author: Lukas Siedentop
 */

#include "nodelist.hpp"

using namespace std;

/* Nodelist */

nodelist::nodelist() :
		min(coordinate()), max(coordinate()), periodic(0), name("") {
}

nodelist::nodelist(bool periodicity, string n) :
		min(coordinate()), max(coordinate()), periodic(periodicity), name(n) {
}

nodelist::nodelist(int pattern, bool periodicity) :
		min(coordinate(-5, -5, -5)), max(coordinate(5, 5, 5)), periodic(
				periodicity) {

	// TODO: periodische RBD

	switch (pattern) {
	case 1: { // Zufallsmuster
		cout << "Generiere Zufallsmuster mit 1000 Punkten in 10^3 Box..."
				<< endl;

		name = "Zufallsmuster";

		// Seed für Zufallsgenerator
		srand(time(NULL));
		for (unsigned n = 0; n < 1000; n++) {
			node* tmp = new node(this,
					(coordinate(3) * 10) - coordinate(5, 5, 5));
			tmp->setEdgenode(1);
			this->push_back(tmp);
		}
		cout << "Zufallsmuster";
		break;
	}
	case 2: { // Diamantmuster (Dirks method)
		cout << "Generiere Diamantmuster mit 1000 Punkten in 10^3 Box..."
				<< endl;

		name = "Diamantmuster";

		int x = -5, y = -5, z = -5;
		int numx = 5, numy = 5, numz = 5;

		while (numz - z > 0) {
			y = -5;
			while (numy - y > 0) {
				x = -5;
				while (numx - x > 0) {
					if (not ((x + y + z) % 2)) {
						// erster Punkt: (0.25,0.25,0.25)
						this->push_back(
								new node(this,
										(coordinate(0.25, 0.25, 0.25)
												+ coordinate(x, y, z))));
						// zweiter Punkt: (0.75,0.75,0.75)
						this->push_back(
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
		for (nodelist::iterator n = this->begin(); n != this->end(); ++n) {
			(*n)->setEdgenode(1);
		}

		cout << "Diamantmuster";
		break;
	}
	}
	// Nachbarn generieren
	setNeighbours();

	cout << " generiert. Statistik:" << endl << listStats() << endl;
}

nodelist::~nodelist() {
	// Periodizität des Musters
	//periodic;

	// Extremalwerte
	min.~coordinate();
	max.~coordinate();

	//TODO nodes selbst mit delete löschen nötig?
	//for (vector<node*>::iterator n = this->begin(); n!=this->end(); ++n) {
	//	(*n)->~node();
	//}

	// Vektor
	vector<node*>().swap(*this);
	//this->~vector();
	//this = vector<node*>();
}

/*
 vector<class node> nodelist::getVector() const {
 vector<class node> nodes;
 for (unsigned i = 0; i < this->size(); i++) {
 nodes.push_back(*this->at(i));
 }
 return nodes;
 }
 */

void nodelist::setNeighbours() {
	// TODO: periodisch
	// vier nächsten Knoten bekommen
	for (unsigned neighs = 0; neighs < 4; neighs++) {
		// Knoteniteration
		for (nodelist::iterator n = begin(); n != end(); ++n) {
			node* neigh = NULL;
			double distSqr = std::numeric_limits<double>::infinity();
			// 2. Knoteniteration
			for (nodelist::iterator nextn = begin(); nextn != end(); ++nextn) {
				// wenn er schon Nachbar ist, nächster
				if ((*n)->isNeighbour((*nextn)) || n == nextn) {
					continue;
				}
				// Differenzvektor
				coordinate diff = *((*nextn)->getPosition())
						- *((*n)->getPosition());
				if (diff.lengthSqr() < distSqr) {
					distSqr = diff.lengthSqr();
					neigh = *nextn;
				}
			}
			if (neigh) {
				(*n)->addNeighbour(neigh);
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
	// Extremalwerte anpassen
	min += shifter;
	max += shifter;

	// Einträge shiften
	for (nodelist::iterator it = begin(); it != end(); ++it) {
		(*it)->shift(shifter);
	}
}

void nodelist::scaleList(double a) {
	// Extremalwerte anpassen
	min *= a;
	max *= a;

	// Einträge shiften
	for (vector<class node*>::iterator it = begin(); it != end(); ++it) {
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
	return size() / getVolume();
}

double nodelist::getVolume() {
	// Volumen berechnen
	double volume = 1;

	coordinate lengths = getLengths();

	for (unsigned int i = 0; i < lengths.dimensions(); i++) {
		volume *= lengths[i];
	}
	return volume;
}

node* nodelist::getAt(coordinate point) {
	// Iteration über alle Knoten
	for (nodelist::iterator it = begin(); it != end(); ++it) {
		if (periodic) {
			coordinate pos = (*(*it)->getPosition());
			//coordinate diff = point - pos;
			if ((std::min(fabs(point[0] - pos[0]),
					(fabs(fabs(point[0] - pos[0]) - getLengths()[0])))
					< 0.000001)
					&& (std::min(fabs(point[1] - pos[1]),
							(fabs(fabs(point[1] - pos[1]) - getLengths()[1])))
							< 0.000001)
					&& (std::min(fabs(point[2] - pos[2]),
							(fabs(fabs(point[2] - pos[2]) - getLengths()[2])))
							< 0.000001)) {
				return (*it);
			}
		}
		if (!periodic && (*(*it)->getPosition() == point)) {
			return (*it);
		}
	}
	return NULL;
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
	for (vector<class node *>::iterator n = begin(); n != end(); ++n) {
		if ((*n)->isEdgenode()) {
			ctr++;
		}
	}
	return ctr;

}

void nodelist::display() {
	std::cout << "Number of elements: " << size() << std::endl;

	int i = 1;

	// Knoteniteration
	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		cout << char(9) << i << "-th element: " << (*(*it)->getPosition())
				<< " with " << (*it)->countNeighbours() << " neighbours at: ";
		// Nachbarniteration
		for (vector<class node *>::iterator neighIt =
				(*it)->getNeighbours()->begin();
				neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
			cout << (*(*neighIt)->getPosition()) << " ";
		}
		cout << endl;
		// weiterzählen
		i++;
	}
}

vector<coordinate> nodelist::getShifters() {
	vector<coordinate> shifters;

	// Boxlängen
	double lx = getLengths()[0], ly = getLengths()[1], lz = getLengths()[2];

	// alle Kombinationen durchgehen
	for (int ix = -1; ix < 2; ix++) {
		for (int iy = -1; iy < 2; iy++) {
			for (int iz = -1; iz < 2; iz++) {
				// die Verschiebung hinzufügen
				shifters.push_back(coordinate(ix * lx, iy * ly, iz * lz));
			}
		}
	}
	return shifters;
}

vector<coordinate> nodelist::getShifted(coordinate mid, double halfExtend) {
	vector<coordinate> shifters;

	// Boxlängen
	double lx = getLengths()[0], ly = getLengths()[1], lz = getLengths()[2]; // 5

	// alle Kombinationen durchgehen, aber nur wenn eine Seitenlänge der gegebenen Box über den Rand geht ...
	for (int ix = -1; ix < 2; ix++) {
		if ((fabs(mid[0] + (ix * halfExtend)) >= fabs(ix * (lx / 2.0)))) {
			for (int iy = -1; iy < 2; iy++) {
				if ((fabs(mid[1] + (iy * halfExtend)) >= fabs(iy * (ly / 2.0)))) {
					for (int iz = -1; iz < 2; iz++) {
						if ((fabs(mid[2] + (iz * halfExtend))
								>= fabs(iz * (lz / 2.0)))) {
							// ... die Verschiebung hinzufügen (Kugel entgegengesetzt zu Muster verschieben!)
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
	// Zähler
	int ctr = 0;

	// Iteration über alle Punkte
	for (vector<class node*>::iterator it = begin(); it != end(); ++it) {
		// Boundingbox
		if ((fabs((*(*it)->getPosition())[0] - mid[0]) < r)
				&& (fabs((*(*it)->getPosition())[1] - mid[1]) < r)
				&& (fabs((*(*it)->getPosition())[2] - mid[2]) < r)) {
			if ((*((*it)->getPosition()) - mid).lengthSqr() < rSqr) {
				ctr++;
			}
		}
	}
	return ctr;
}

int nodelist::pointsInsidePeriodic(coordinate mid, double r) {
	// Zähler
	int ctr = 0;

	// Verschiebevektoren
	vector<coordinate> shifted = getShifted(mid, r);

	// Abstandsvergleich mit in jeder Raumrichtung verschobenen Kugel
	for (unsigned int s = 0; s < shifted.size(); s++) {
		// Iteration über alle Punkte
		for (vector<class node*>::iterator it = begin(); it != end(); ++it) {
			// Boundingbox
			if ((fabs((*(*it)->getPosition())[0] + shifted[s][0]) < r)
					&& (fabs((*(*it)->getPosition())[1] + shifted[s][1]) < r)
					&& (fabs((*(*it)->getPosition())[2] + shifted[s][2]) < r)) {

				// Würfel sind leider nicht strikt konvex... (Local density fluctuations, hyperuniformity, and order metrics. Salvatore Torquato et. al.)
				//ctr++;

				// liegt der Punkt in der Kugel?
				if ((((*(*it)->getPosition()) + shifted[s]).lengthSqr())
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
	stream << commentDelimeter << "Number of nodes: " << size() << endl;
	stream << commentDelimeter << "Bounding box: " << min << ", " << max
			<< endl;
	stream << commentDelimeter << "Mid of box: " << getMid() << endl;
	stream << commentDelimeter << "Volume: " << getVolume() << endl;
	stream << commentDelimeter << "Density of points: " << getDensity() << endl;

	return stream.str();
}

double nodelist::normalize() {
	// Dichte auf 1 setzen
	double factor = pow(getDensity(), 1. / 3.);
	//double factor = ((M_PI * pow(300, 2./3.)) / 6) / getMaxFeatureSize();//pow(getDensity(), 1. / 3.);
	//scaleList(factor);
	setDensity(1);

	// Mittelpunkt anpassen
	shiftList(getMid() * -1);
	/*
	 coordinate* shifter = (*this)[0]->getPosition();
	 // Dem Ursprung nächsten Knoten finden
	 for (vector<class node*>::iterator it = begin(); it != end(); ++it) {
	 if (shifter->lengthSqr() > (*it)->getPosition()->lengthSqr()) {
	 shifter = (*it)->getPosition();
	 }
	 }
	 shiftList(*shifter * -1);
	 */
	return factor;
}

/**
 * Schreibt die Anzahl der Nachbarn jedes Knotens untereinander in eine Datei.
 */
vector<double> nodelist::neighbourDistribution() {
	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iterationsvariablen
	int counter;

	// Iteration über jeden Knoten
	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		if (!(*it)->isEdgenode()) {
			counter = 0;
			for (vector<class node *>::iterator neighIt =
					(*it)->getNeighbours()->begin();
					neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
				// Nachbarn zählen
				counter++;
			}

			// Wert schreiben
			data.push_back(counter);
		}
	}

	return data;
}

/**
 * Schreibt alle Distanzen zwischen den Knoten und ihrer Nachbarn in eine Datei.
 */
vector<double> nodelist::lengthDistribution(bool plot) {
	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iteration über jeden Knoten
	for (vector<class node *>::iterator n = begin(); n != end(); ++n) {
		if (!(*n)->isEdgenode()) {
			// über Nachbarn iterieren
			for (vector<class node *>::iterator nn =
					(*n)->getNeighbours()->begin();
					nn != (*n)->getNeighbours()->end(); ++nn) {
				// je nach Periodizität Wert schreiben
				if (periodic) {
					data.push_back((*n)->euklidianPeriodic((*nn)));
				} else {
					data.push_back((*n)->euklidian((*nn)));
				}
			}
		}
	}

	// Daten sortieren...
	data = mergeSort(data);
	vector<double> halfData;

	// ...und jeden zweiten behalten.
	for (unsigned int i = 0; i < data.size(); i += 2) {
		halfData.push_back(data[i]);
	}

	return halfData;
}

/**
 * Schreibt die Winkel zwischen den Nachbarn in eine Datei.
 */
vector<double> nodelist::angleDistribution() {
	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iteration über jeden Knoten
	for (vector<node*>::iterator nodeIter = begin(); nodeIter != end();
			++nodeIter) {
		if (!(*nodeIter)->isEdgenode()) {
			// über Nachbarn iterieren
			for (vector<node*>::iterator neighIter1 =
					(*nodeIter)->getNeighbours()->begin();
					neighIter1 != (*nodeIter)->getNeighbours()->end();
					++neighIter1) {

				// über Nachbarn des Nachbarn iterieren
				for (vector<node*>::iterator neighIter2 =
						(*neighIter1)->getNeighbours()->begin();
						neighIter2 != (*neighIter1)->getNeighbours()->end();
						++neighIter2) {
					// wenn nicht der Winkel zwischen sich berechnet werden soll...
					if (!(*nodeIter)->equals(*neighIter2)) {
						// ...Wert berechnen & schreiben, umgerechnet in Grad
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

	// Daten sortieren...
	data = mergeSort(data);
	vector<double> halfData;

	// ...und jeden zweiten behalten.
	for (unsigned int i = 0; i < data.size(); i += 2) {
		halfData.push_back(data[i]);
	}

	return halfData;
}

/**
 * Calculates the hyperuniformity of the pattern saved in the list.
 */
void nodelist::hyperuniformity(vector<vector<double> >& variance) {
	// number of radii
	unsigned int nr = 50;
	// radiusincrement + maximal radius
	double rMax = periodic ? getLengths().min() : getLengths().min() / 3, dr =
			rMax / nr;

	// number of spheres/iterations
	unsigned int n = 500;

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
	cout << "Hat " << t << " Clicks gedauert ("
			<< (((float) t) / CLOCKS_PER_SEC) << "s)" << endl;

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

	//return variance;
}

/**
 * Schreibt eine pov-datei der Stäbe
 */
void nodelist::writePOV() {

	const char outfileName[] = "./data/staebe.pov";

	vector<string> data;

	for (nodelist::iterator nodeIter = begin(); nodeIter != end(); ++nodeIter) {
		stringstream stream;

		// Kugel am Gelenk
		stream << "sphere{"
				<< (*nodeIter)->getPosition()->toString("<", ",", ">") << ",r}"
				<< endl;

		for (vector<node*>::iterator neighIter =
				(*nodeIter)->getNeighbours()->begin();
				neighIter != (*nodeIter)->getNeighbours()->end(); ++neighIter) {

			//  nur wenn nachbar und knoten nicht über die grenze verbunden sind...
			if (isPeriodic()
					&& (*nodeIter)->euklidian((*neighIter))
							< getMaxFeatureSize()) {
				// Zylinder von A nach B TODO: den zurück nicht...
				stream << "cylinder{"
						<< (*nodeIter)->getPosition()->toString("<", ",", ">")
						<< ","
						<< (*neighIter)->getPosition()->toString("<", ",", ">")
						<< ",r}" << endl;

			} else if (!isPeriodic()) {
				// Zylinder von A nach B TODO: den zurück nicht...
				stream << "cylinder{"
						<< (*nodeIter)->getPosition()->toString("<", ",", ">")
						<< ","
						<< (*neighIter)->getPosition()->toString("<", ",", ">")
						<< ",r}" << endl;
			}
		}
		data.push_back(stream.str());
	}

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Daten schreiben
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i];
	}

	cout << "Stäbe in " << outfileName << " geschrieben." << endl;

}
