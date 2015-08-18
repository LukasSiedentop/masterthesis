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
	double factor = density / getDensity();
	scaleList(factor);
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
	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		// Vergleich der Koordinaten
		if (*(*it)->getPosition() == point) {
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

double nodelist::getMaxFeatureSize() {
	return getLengths().length() / 2;
}

void nodelist::display() {
	std::cout << "Anzahl Elemente: " << size() << std::endl;

	int i = 1;

	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {

		cout << char(9) << i << "-tes Element: " << (*it)->getPosition()
				<< " mit " << (*it)->countNeighbours() << " Nachbarn bei: ";

		for (vector<class node *>::iterator neighIt =
				(*it)->getNeighbours()->begin();
				neighIt != (*it)->getNeighbours()->end(); ++neighIt) {
			cout << (*neighIt)->getPosition();
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
							// ... die Verschiebung hinzufügen
							shifters.push_back(
									mid
											+ coordinate(ix * lx, iy * ly,
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

// TODO: paralellisieren TODO: prüfen
int nodelist::pointsInside(coordinate mid, double r) {
	// Zähler
	int ctr = 0;

	// Iteration über alle Punkte
	for (vector<class node *>::iterator it = begin(); it != end(); ++it) {
		if ((*((*it)->getPosition()) - mid).lengthSqr() < r * r) {
			ctr++;
		}
	}
	return ctr;
}

// TODO: paralellisieren TODO: noch falsch
int nodelist::pointsInsidePeriodic(coordinate mid, double r) {
	// Zähler
	int ctr = 0;

	// Verschiebevektoren
	//vector<coordinate> shifters = getShifters(); // funktioniert
	vector<coordinate> shifted = getShifted(mid, r); // funktioniert nicht

	// Abstandsvergleich mit in jeder Raumrichtung verschobenen Kugel
	for (unsigned int s = 0; s < shifted.size(); s++) {

		// Iteration über alle Punkte
		for (vector<class node *>::iterator it = begin(); it != end(); ++it) {

			/*
			 // Box statt Kugel
			 if ((( tmp + shifters[s] )[0] < r)&&( tmp + shifters[s] )[1] > r)&&(( tmp + shifters[s] )[2] > r)) {
			 ctr++;
			 }
			 }
			 */

			//if (((*((*it)->getPosition()) - mid) - shifters[s]).length() < r) {
			//if ((*it)->euklidian((mid + shifters[s])) < r) {
			//cout << "r: " << r <<" p: " << (*(*it)->getPosition()) << " s: " << shifted[s] << " p->euklidian(s): " << (*it)->euklidian(shifted[s]) << endl;
			if ((((*(*it)->getPosition()) + shifted[s]).lengthSqr())
					< (r * r)) {
				ctr++;
				//cout << "r: " << r << " r*r: " << (r*r) <<" p: " << (*(*it)->getPosition()) << " s: " << shifted[s] << " p->euklidian(s)^2: " << (((*(*it)->getPosition()) + shifted[s]).lengthSqr()) << endl;
			}
		}
	}
	return ctr;
}

string nodelist::listStats(const char commentDelimeter[]) {
	stringstream stream;
	stream << commentDelimeter << "Anzahl Knoten: " << size() << endl;
	stream << commentDelimeter << "Extremalwerte: " << min << ", " << max
			<< endl;
	stream << commentDelimeter << "Volumen: " << getVolume() << endl;
	stream << commentDelimeter << "Punktdichte: " << getDensity() << endl;

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
	plotHist(halfData, 0.6, 1.05, 25, "Länge");
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
	const int n = 20;

	// Fortschrittsbalken
	boost::progress_display show_progress(n);

	// Datenstruktur um Anzahl der Punkte in Kugel zu speichern
	vector<vector<double> > data;

	unsigned const int nr = ceil(rMax / dr);

	data.resize(nr + 1);
	for (unsigned int i = 0; i < nr; ++i) {
		data[i].resize(n);
	}

	// erstelle n Kugeln...
	for (int i = 0; i < n; i++) {
		// Zufälligen Mittelpunkt auswählen (ist ~2 clicks schneller als einzeln die koordinaten generieren)
		// TODO: Rand beachten bei nicht periodisch -> zufälligen Mittelpunkt wählen, sodass kugel immer im muster liegt
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

	// Varianz über Radius plotten
	plot2D(variance, "Radius R", "Varianz  {/Symbol s}^2(R)");

}
