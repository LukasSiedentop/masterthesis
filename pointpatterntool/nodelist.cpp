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
		min(NULL), max(NULL), periodic(0) {
}

nodelist::nodelist(bool periodicity) :
		min(NULL), max(NULL), periodic(periodicity) {
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
	for (vector<class node*>::iterator it = begin(); it != end(); ++it) {
		node* node = it;
		node->shift(shifter);
	}
}

void nodelist::scaleList(double a, double b, double c) {
	// Extremalwerte anpassen
	min *= coordinate(a, b, c);
	max *= coordinate(a, b, c);

	// Einträge shiften
	for (vector<class node*>::iterator it = begin(); it != end(); ++it) {
		node* node = it;
		node->scale(a, b, c);
	}
}

bool nodelist::isPeriodic() {
	return periodic;
}

void nodelist::setDensity(double density) {
	double factor = density / getDensity();
	scaleList(factor, factor, factor);
}

double nodelist::getDensity() {
	return size() / getVolume();
}

double nodelist::getVolume() {
	// Volumen berechnen
	double volume = 1;
	for (vector<double>::iterator it = getLengths().begin();
			it != getLengths().end(); ++it) {
		volume *= it;
	}
	return volume;
}

node * nodelist::getAt(coordinate point) {
	// Iteration über alle Knoten
	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		node* node = it;

		// Punkt in Box holen
		/*
		 if (periodic) {
		 node = coordinate::getVec(coordinate(0,0,0), node->getPosition(), getShifters());
		 }
		 */

		// Vergleich der Koordinaten
		if (node->getPosition() == point) {
			return node;
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

void nodelist::display() {
	std::cout << "Anzahl Elemente: " << size() << std::endl;

	int i = 1;

	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		node* n = it;

		cout << char(9) << i << "-tes Element: " << n->getPosition() << " mit "
				<< n->countNeighbours() << " Nachbarn bei: ";

		for (vector<class node *>::iterator itNeigh =
				n->getNeighbours()->begin();
				itNeigh != n->getNeighbours()->end(); ++itNeigh) {
			node* neigh = itNeigh;
			cout << neigh->getPosition();
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
				coordinate coord;
				coord.push_back(ix * lx);
				coord.push_back(iy * ly);
				coord.push_back(iz * lz);

				shifters.push_back(coord);
			}
		}
	}

	return shifters;
}

int nodelist::pointsInside(double r, coordinate mid) {

	int ctr = 0;
	vector<coordinate> shifters = getShifters();
	// TODO: paralelisieren
	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		node* n = it;

		if (periodic) {
			// Abstandsvergleich mit in jeder Raumrichtung verschobenen Kugel
			for (unsigned int i = 0; i < shifters.size(); i++) {
				if (n->euklidian(mid + shifters[i]) < r) {
					ctr++;
				}
			}
		} else if (n->euklidian(mid) < r) {
			ctr++;
		}
	}

	return ctr;
}

string nodelist::listStats(const char commentDelimeter[]) {
	double volume = getVolume();
	double density = getDensity();

	stringstream stream;
	stream << commentDelimeter << "Anzahl Knoten: " << size << endl;
	stream << commentDelimeter << "Extremalwerte: " << min << ", " << max
			<< endl;
	stream << commentDelimeter << "Volumen: " << volume << endl;
	stream << commentDelimeter << "Punktdichte: " << density << endl;

	return stream.str();
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
		node* n = it;
		counter = 0;
		for (vector<class node *>::iterator neighIt = n->getNeighbours()->begin(); neighIt != n->getNeighbours()->end(); ++neighIt) {
				node* neigh = neighIt;
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
