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
		min(coordinate()), max(coordinate()), periodic(0) {
}

nodelist::nodelist(bool periodicity) :
		min(coordinate()), max(coordinate()), periodic(periodicity) {
}

nodelist::nodelist(vector<node>& vec, bool periodicity) :
		periodic(periodicity) {
	// Extremalwerte bestimmen
	double inf = numeric_limits<double>::infinity();
	double minX = inf, minY = inf, minZ = inf, maxX = -inf, maxY = -inf, maxZ =
			-inf;

	//coordinate minimum({inf, inf, inf},3);
	//coordinate maximum({-inf, -inf, -inf},3);

	for (unsigned n = 0; n < vec.size(); n++) {

		minX = std::min(minX, (*vec[n].getPosition())[0]);
		minY = std::min(minY, (*vec[n].getPosition())[1]);
		minZ = std::min(minZ, (*vec[n].getPosition())[2]);

		maxX = std::max(maxX, (*vec[n].getPosition())[0]);
		maxY = std::max(maxY, (*vec[n].getPosition())[1]);
		maxZ = std::max(maxZ, (*vec[n].getPosition())[2]);

		// Liste des Knotens neu setzen
		//vec[n].setList(this);
		// Knoten dieser Liste hinzufügen
		this->push_back(new node(vec[n], this));
	}

	setMins(coordinate(minX, minY, minZ));
	setMaxs(coordinate(maxX, maxY, maxZ));
}

nodelist::~nodelist() {
	vector<node*>().swap(*this);

	min.~coordinate();
	max.~coordinate();

	// Periodizität des Musters
	//periodic;
}

vector<class node> nodelist::getVector() {
	vector<class node> nodes;
	for (unsigned i = 0; i < this->size(); i++) {
		nodes.push_back(*this->at(i));
	}
	return nodes;
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

bool nodelist::isPeriodic() {
	return periodic;
}

void nodelist::setDensity(double density) {
	scaleList(pow((getDensity() / density), 1. / 3.));
}

double nodelist::getDensity() {
	return size() / getVolume();
}

double nodelist::getVolume() {
	// Volumen berechnen
	double volume = 1;

	coordinate lengths = getLengths();

	for (coordinate::iterator it = lengths.begin(); it != lengths.end(); ++it) {

		volume *= (*it);
	}
	return volume;
}

node * nodelist::getAt(coordinate point) {
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

		// Vergleich der Koordinaten
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

void nodelist::display() {
	std::cout << "Anzahl Elemente: " << size() << std::endl;

	int i = 1;

	// Knoteniteration
	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		cout << char(9) << i << "-tes Element: " << (*(*it)->getPosition())
				<< " mit " << (*it)->countNeighbours() << " Nachbarn bei: ";
		// Nachbarniteration
		for (vector<class node *>::iterator neighIt =
				(*it)->getNeighbours()->begin();
				neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
			cout << (*(*neighIt)->getPosition());
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
		if ((fabs(mid[0] + (ix * halfExtend)) > fabs(ix * (lx / 2.0)))) {
			for (int iy = -1; iy < 2; iy++) {
				if ((fabs(mid[1] + (iy * halfExtend)) > fabs(iy * (ly / 2.0)))) {
					for (int iz = -1; iz < 2; iz++) {
						if ((fabs(mid[2] + (iz * halfExtend))
								> fabs(iz * (lz / 2.0)))) {
							// ... die Verschiebung hinzufügen (Kugel entgegengesetzt zu Muster verschieben!)
							shifters.push_back(
									mid
											- coordinate(ix * lx, iy * ly,
													iz * lz));
							//cout << "r: " << halfExtend << " s: " << mid+coordinate(ix * halfExtend, iy * halfExtend, iz * halfExtend) << " p: " << coordinate(ix * lx, iy * ly, iz * lz) << endl;
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
	stream << commentDelimeter << (periodic ? "Periodisch" : "Nicht periodisch")
			<< endl;
	stream << commentDelimeter << "Anzahl Knoten: " << size() << endl;
	stream << commentDelimeter << "Extremalwerte: " << min << ", " << max
			<< endl;
	stream << commentDelimeter << "Mitte: " << getMid() << endl;
	stream << commentDelimeter << "Volumen: " << getVolume() << endl;
	stream << commentDelimeter << "Punktdichte: " << getDensity() << endl;

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
void nodelist::neighbourDistribution() {
	const char outfileName[] = "./data/statistics/neighbourDistribution.dat";
	cout << "Zähle die Nachbarn jedes Knotens..." << endl;

	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iterationsvariablen
	int counter;

	// Iteration über jeden Knoten
	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		counter = 0;
		for (vector<class node *>::iterator neighIt =
				(*it)->getNeighbours()->begin();
				neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
			// Nachbarn zählen
			// TODO: wenn nicht periodisch nur nachbarn die nicht am Rand sind.
			counter++;
		}

		// Wert schreiben
		data.push_back(counter);
	}

	// Daten schreiben
	writeHist(data, true, "# Anzahl Nächster Nachbarn", outfileName);

	// Statistiken
	cout << "Nachbarstatistik:" << endl;
	cout << stats(data);

	// Daten ploten
	plotHist(data, 0, 10, 10, "Anzahl Nächster Nachbarn");
}

/**
 * Schreibt alle Distanzen zwischen den Knoten und ihrer Nachbarn in eine Datei.
 */
void nodelist::lengthDistribution() {
	const char outfileName[] = "./data/statistics/lenghtDistribution.dat";
	cout << "Bestimme die Längen zwischen benachbarten Punkten..." << endl;

	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iteration über jeden Knoten
	for (vector<class node *>::iterator n = begin(); n != end(); ++n) {
		// über Nachbarn iterieren
		for (vector<class node *>::iterator nn = (*n)->getNeighbours()->begin();
				nn != (*n)->getNeighbours()->end(); ++nn) {
			// je nach Periodizität Wert schreiben
			if (periodic) {
				data.push_back((*n)->euklidianPeriodic((*nn)));
			} else {
				data.push_back((*n)->euklidian((*nn)));
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

	// Daten schreiben
	writeHist(halfData, true, "# Längen", outfileName);

	// Statistiken
	cout << "Längenstatistik:" << endl;
	cout << stats(halfData);

	// Daten ploten
	plotHist(halfData, floor((*halfData.begin())) - 0.5,
			ceil((*halfData.end())) + 0.5, 30, "Länge");

}

/**
 * Schreibt die Winkel zwischen den Nachbarn in eine Datei.
 */
void nodelist::angleDistribution() {
	const char outfileName[] = "./data/statistics/angleDistribution.dat";

	cout << "Bestimme die Winkel zwischen benachbarten Punkten..." << endl;

	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iteration über jeden Knoten
	for (vector<node*>::iterator nodeIter = begin(); nodeIter != end();
			++nodeIter) {
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
												(*nodeIter), (*neighIter2)));
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

	// Daten sortieren...
	data = mergeSort(data);
	vector<double> halfData;

	// ...und jeden zweiten behalten.
	for (unsigned int i = 0; i < data.size(); i += 2) {
		halfData.push_back(data[i]);
	}

	// Daten schreiben
	writeHist(halfData, true, "# Winkel", outfileName);

	// Statistiken
	cout << "Winkelstatistik:" << endl;
	cout << stats(halfData);

	// Daten ploten
	plotHist(halfData, 0, 180, 180, "Winkel");

	cout << "Winkeldaten in " << outfileName << " geschrieben." << endl;

}

/**
 * Errechnet die Hyperuniformität des Musters.
 */
void nodelist::hyperuniformity() {

	const char outfileName[] = "./data/statistics/hyperuniformity.dat";

	// Radiusincrement
	double dr = 0.2, rMax = 11;

	// Anzahl Kugeln bzw. iterationen
	const int n = 1000;

	// Datenstruktur um Anzahl der Punkte in Kugel zu speichern
	vector<vector<double> > data;

	unsigned int nr = ceil(rMax / dr);

	data.resize(nr + 1);
	for (unsigned int i = 0; i < nr; ++i) {
		data[i].resize(n);
	}

	// Seed für Zufallsgenerator
	srand(time(NULL));
	if (periodic) {

		// Fortschrittsbalken
		boost::progress_display show_progress(n);

		// erstelle n Kugeln...
		for (int i = 0; i < n; i++) {
			// Zufälligen Mittelpunkt auswählen
			coordinate mid(3);
			mid *= getLengths();
			mid -= getLengths() / 2;

			// Über Radius iterieren
			for (unsigned int j = 0; j < nr; j++) {
				data[j][i] = pointsInsidePeriodic(mid, j * dr);
			}

			// Fortschritt
			++show_progress;
		}
	} else {

		// Fortschrittsbalken
		boost::progress_display show_progress(nr);

		double r;
		double rSqr;
		// Über Radius iterieren
		for (unsigned int j = 0; j < nr; j++) {
			r = j * dr;
			rSqr = r * r;
			for (int i = 0; i < n; i++) {
				// zufälligen Mittelpunkt wählen, sodass Kugel immer im Muster liegt
				coordinate mid(3);
				coordinate range = getLengths()
						- coordinate(2 * r, 2 * r, 2 * r);
				mid *= range;
				mid -= range / 2;

				data[j][i] = pointsInside(mid, r, rSqr);
			}

			// Fortschritt
			++show_progress;
		}
	}

	//cout << "Radius" << char(9) << "Volumen" << char(9) << "Erwartungswert" << char(9) << "Verhältnis" << endl;
	// Erwartungswert aller Radii berechnen
	// Der sollte genau die Dichte des Musters mal das Kugelvolumen sein!
	vector<double> expectedValue;
	expectedValue.resize(nr + 1);
	// Radiusiteration
	for (unsigned int j = 0; j < nr; j++) {
		// Kugeliteration
		for (int i = 0; i < n; i++) {
			expectedValue[j] += data[j][i];
		}
		expectedValue[j] = expectedValue[j] / n;
		//cout << dr*j << char(9) << 4/3 * M_PI * pow(dr*j, 3) << char(9) << expectedValue[j] << char(9) << expectedValue[j]/(4/3 * M_PI * pow(dr*j, 3)) << endl;
	}

	// Varianz für jeden Radius berechnen Gnuplot brauchts so: variance[i][j], (i,j) => (rows,colums)
	vector<vector<double> > variance;
	variance.resize(nr + 1);
	// Radiusiteration
	for (unsigned int j = 0; j < nr; j++) {
		variance[j].resize(2);
		variance[j][0] = j * dr;

		// Kugeliteration
		for (int i = 0; i < n; i++) {
			variance[j][1] += pow((data[j][i] - expectedValue[j]), 2);
		}
		variance[j][1] = variance[j][1] / n;
	}

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Daten schreiben
	outfile << "# Radius" << char(9) << "Varianz" << endl;
	for (unsigned int j = 0; j < nr; j++) {
		outfile << variance[j][0] << char(9) << variance[j][1] << endl;
	}

	double fitMax = std::min(getLengths()[0],
			std::min(getLengths()[1], getLengths()[2]));

	// Varianz über Radius plotten
	plot2D(variance, 0, dr, rMax, periodic ? : fitMax / 2, "Radius R",
			"Varianz  {/Symbol s}^2(R)");
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
