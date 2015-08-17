/*
 * pointpatterntool.cpp
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <ctime>

#include <voro++/voro++.hh>

#include <boost/tuple/tuple.hpp>
#include <boost/progress.hpp>

#include "functions.hpp"
#include "nodelist.hpp"

#include "coordinate.hpp"

using namespace std;

/**
 * Schreibt die Optionen auf die Koinsole und wartet auf einen Input.
 */
int gui() {
	cout << "---------------------------------------------------" << endl;
	cout << "1 - Vert. der Anzahl der Nachbarn" << endl;
	cout << "2 - Vert. der Abstände zu den Nachbarn" << endl;
	cout << "3 - Vert. der Winkel der Nachbarn" << endl;
	cout << "4 - Grad der Hyperuniformity" << endl;

	cout << "5 - Liste darstellen" << endl;
	cout << "6 - Listenstatistik" << endl;
	cout << "7 - Muster gnuplotten" << endl;
	cout << "8 - POV-Ray Datei schreiben" << endl;

	cout << "9 - Muster schreiben (voro++)" << endl;
	cout << "10 - Voro++ ausprobieren" << endl;

	cout << "11 - weitere Liste einlesen" << endl;
	cout << "12 - Listen vergleichen" << endl;

	cout << "0 - Nichts." << endl;
	cout << "---------------------------------------------------" << endl;
	cout << "Was möchtest du über die Punkte wissen? (0-12; default: 0) << ";

	int option = 0;

	return input(option);
}

/**
 * Liest zwei Dateien ein. Die erste enthält Punkte im R^3, die zweite Nachbarn dieser Punkte. Beispiel:
 * nodeFile:	neighbourFile:
 * 1 2 3		2 3 4
 * 1 2 3		3 4 5
 * 2 3 4		1 2 3
 * 2 3 4		3 4 5
 * 3 4 5		1 2 3
 * 3 4 5		2 3 4
 */
nodelist * readfile(const char * nodes, const char * neighbours) {
	cout << "Lese Knotendatei " << nodes << " und Nachbardatei " << neighbours
			<< " ein." << endl;

	// input filestream beider Dateien
	ifstream nodeFile, neighbourFile;
	nodeFile.open(nodes, std::ifstream::in);
	neighbourFile.open(neighbours, std::ifstream::in);

	cout << "Ist das Muster periodisch? (0,1; default: 0) << ";

	bool periodic = input(0);

	// Listenkopf
	nodelist * list = new nodelist(periodic);

	if (periodic) {
		cout << "Nehme eine Box von (-5,5)^3 an..." << endl;
		list->setMins(coordinate(-5, -5, -5));
		list->setMaxs(coordinate(5, 5, 5));
	}

	// Position des Punktes
	double x, y, z;
	// Position des Nachbarn
	double nx, ny, nz;

	// Extremalwerte
	double minX = 0, minY = 0, minZ = 0, maxX = 0, maxY = 0, maxZ = 0;

	// einlesen der Dateien
	while ((nodeFile >> x >> y >> z) && (neighbourFile >> nx >> ny >> nz)) {

		node * n = list->getAt(coordinate(x, y, z));
		// Wenn an der Knoten Position noch kein Knoten existiert...
		if (!n) {
			// ... füge ihn der Liste hinzu
			n = new node(list, x, y, z);
			list->push_back(n);
		}

		node * neighbour = list->getAt(coordinate(nx, ny, nz));
		// Wenn an der Nachbarknoten Position noch kein Knoten existiert...
		if (!neighbour) {
			// ... füge ihn der Liste hinzu
			neighbour = new node(list, nx, ny, nz);
			list->push_back(neighbour);
		}

		// die Nachbarschaft setzen
		n->getNeighbours()->push_back(neighbour);
		neighbour->getNeighbours()->push_back(n);

		// wenn nicht periodisch...
		if (!periodic) {
			// ... Extremalwerte bekommen
			minX = min(minX, x);
			minY = min(minY, y);
			minZ = min(minZ, z);

			maxX = max(maxX, x);
			maxY = max(maxY, y);
			maxZ = max(maxZ, z);
		}
	}

	// Extremalwerte setzen
	if (!periodic) {
		list->setMins(coordinate(minX, minY, minZ));
		list->setMaxs(coordinate(maxX, maxY, maxZ));
	}

	cout << "Dateien eingelesen und Liste mit " << list->size()
			<< " Knoten erstellt." << endl;
	return list;
}/*

nodelist * readfileGui() {
	string nodes = "data/Left.dat", neighbours = "data/Right.dat";

	cout << "Lese zweite Liste ein..." << endl;
	cout << "Knotendatei (default: " << nodes << ") << " << endl;
	nodes = input(nodes);
	cout << "Nachbarndatei (default: " << neighbours << ") << " << endl;
	neighbours = input(neighbours);

	return readfile(nodes.c_str(), neighbours.c_str());
}*/

// TODO: Methoden die nodelist * list brauchen in nodelist?!

/**
 * Schreibt alle Distanzen zwischen den Knoten und ihrer Nachbarn in eine Datei.
 *//*
void lengthDistribution(nodelist * list, const char outfileName[]) {
	cout << "Bestimme die Längen zwischen benachbarten Punkten..." << endl;

	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iterationsvariablen
	Node * nodeIter = list->getFirst();
	Neighbour * neighIter;

	// Iteration über jeden Knoten
	while (nodeIter) {
		neighIter = nodeIter->getNeighbours();
		// über Nachbarn iterieren
		while (neighIter) {
			// je nach Periodizität Wert schreiben
			if (list->isPeriodic()) {
				data.push_back(
						nodeIter->euklidianPeriodic(neighIter->getNode()));
			} else {
				data.push_back(nodeIter->euklidian(neighIter->getNode()));
			}

			// weiter gehts
			neighIter = neighIter->getNext();
		}
		nodeIter = nodeIter->getNext();
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
}*/

/**
 * Schreibt die Winkel zwischen den Nachbarn in eine Datei.
 *//*
void angleDistribution(nodelist * list, const char outfileName[]) {
	cout << "Bestimme die Winkel zwischen benachbarten Punkten..." << endl;

	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iterationsvariablen
	Node *nodeIter = list->getFirst();
	Neighbour *neighIter1, *neighIter2;

	// Iteration über jeden Knoten
	while (nodeIter) {
		neighIter1 = nodeIter->getNeighbours();
		// über Nachbarn iterieren
		while (neighIter1) {
			neighIter2 = neighIter1->getNode()->getNeighbours();
			// über Nachbarn des Nachbarn iterieren
			while (neighIter2) {
				// wenn nicht der Winkel zwischen sich berechnet werden soll...
				if (!nodeIter->equals(neighIter2->getNode())) {
					// ...Wert berechnen & schreiben, umgerechnet in Grad
					if (list->isPeriodic()) {
						data.push_back(
								((180.0) / M_PI)
										* neighIter1->getNode()->anglePeriodic(
												nodeIter,
												neighIter2->getNode()));
					} else {

						data.push_back(
								((180.0) / M_PI)
										* neighIter1->getNode()->angle(nodeIter,
												neighIter2->getNode()));
					}
				}
				// weiter gehts
				neighIter2 = neighIter2->getNext();
			}
			// weiter gehts
			neighIter1 = neighIter1->getNext();
		}
		// weiter gehts
		nodeIter = nodeIter->getNext();
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
}*/

/**
 * Errechnet die Hyperuniformität des Musters.
 *//*
void hyperuniformity(nodelist * list, const char outfileName[]) {

	// Radiusincrement
	double dr = 0.2, rMax = 11;

	// Anzahl Kugeln bzw. iterationen
	const int n = 1000;

	// Fortschrittsbalken
	boost::progress_display show_progress(n);

	// Datenstruktur um Anzahl der Punkte in Kugel zu speichern
	vector<vector<double> > data;
	data.resize(((int) floor(rMax / dr) + 1));
	for (int i = 0; i < ((int) floor(rMax / dr) + 1); ++i) {
		data[i].resize(n);
	}

	for (int i = 0; i < n; i++) {
		// Zufälligen Mittelpunkt auswählen
		double mx = ((double) rand() / RAND_MAX) * (list->lengthX())
				+ list->getMinX();
		double my = ((double) rand() / RAND_MAX) * (list->lengthY())
				+ list->getMinY();
		double mz = ((double) rand() / RAND_MAX) * (list->lengthZ())
				+ list->getMinZ();

		// Über Radius iterieren
		for (int j = 0; j < floor(rMax / dr); j++) {
			// TODO: Rand beachten bei nicht periodisch -> zufälligen Mittelpunkt wählen, sodass kugel immer im muster liegt
			data[j][i] = list->pointsInside(j * dr, mx, my, mz);
		}
		++show_progress;
	}

	//cout << "Radius" << char(9) << "Volumen" << char(9) << "Erwartungswert" << char(9) << "Verhältnis" << endl;
	// Erwartungswert aller Radii berechnen
	// Der sollte genau die Dichte des Musters mal das Kugelvolumen sein!
	vector<double> expectedValue;
	expectedValue.resize(floor(rMax / dr) + 1);
	for (int j = 0; j < floor(rMax / dr); j++) {
		for (int i = 0; i < n; i++) {
			expectedValue[j] += data[j][i];
		}
		expectedValue[j] = expectedValue[j] / n;

		//cout << dr*j << char(9) << 4/3 * M_PI * pow(dr*j, 3) << char(9) << expectedValue[j] << char(9) << expectedValue[j]/(4/3 * M_PI * pow(dr*j, 3)) << endl;

	}

	// Varianz für jeden Radius berechnen Gnuplot brauchts so: variance[i][j], (i,j) => (rows,colums)
	vector<vector<double> > variance;
	variance.resize(floor(rMax / dr) + 1);

	for (int i = 0; i < floor(rMax / dr); i++) {
		variance[i].resize(2);

		variance[i][0] = i * dr;

		for (int j = 0; j < n; j++) {
			variance[i][1] += pow((data[i][j] - expectedValue[i]), 2);
		}
		variance[i][1] = variance[i][1] / n;
	}

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Daten schreiben
	outfile << "# Radius" << char(9) << "Varianz" << endl;
	for (int i = 0; i < floor(rMax / dr); i++) {
		outfile << variance[i][0] << char(9) << variance[i][1] << endl;
	}

	// Varianz über Radius plotten
	plot2D(variance, "Radius R", "Varianz  {/Symbol s}^2(R)");

}*/

/**
 * Bereitet die Daten auf, das Gnuplot die Stäbe plotten kann.
 *//*
void gnuplotPattern(nodelist * list) {
	vector<vector<double> > data;

	Node * nodeIter = list->getFirst();
	while (nodeIter) {
		Neighbour * neighIter = nodeIter->getNeighbours();
		while (neighIter) {
			vector<double> line;
			line.resize(3);

			//  nur wenn nachbar und knoten nicht über die grenze verbunden sind...
			if (list->isPeriodic()
					&& nodeIter->euklidian(neighIter->getNode())
							< list->lengthX() / 2) {
				line[0] = nodeIter->getX();
				line[1] = nodeIter->getY();
				line[2] = nodeIter->getZ();

				data.push_back(line);

				line[0] = neighIter->getNode()->getX();
				line[1] = neighIter->getNode()->getY();
				line[2] = neighIter->getNode()->getZ();

				data.push_back(line);

				line[0] = nan("");
				line[1] = nan("");
				line[2] = nan("");

				data.push_back(line);

				line[0] = nan("");
				line[1] = nan("");
				line[2] = nan("");

				data.push_back(line);
			} else if (!list->isPeriodic()) {
				line[0] = nodeIter->getX();
				line[1] = nodeIter->getY();
				line[2] = nodeIter->getZ();

				data.push_back(line);

				line[0] = neighIter->getNode()->getX();
				line[1] = neighIter->getNode()->getY();
				line[2] = neighIter->getNode()->getZ();

				data.push_back(line);

				line[0] = nan("");
				line[1] = nan("");
				line[2] = nan("");

				data.push_back(line);

				line[0] = nan("");
				line[1] = nan("");
				line[2] = nan("");

				data.push_back(line);
			}

			neighIter = neighIter->getNext();
		}
		nodeIter = nodeIter->getNext();
	}

	// Plotten
	plot3D(data, list->getMinX(), list->getMaxX());
}
*/
/**
 * Schreibt eine pov-datei der Stäbe
 *//*
void writePOV(nodelist * list, const char outfileName[]) {
	vector<string> data;

	node * nodeIter = list->getFirst();
	while (nodeIter) {

		stringstream stream;

		// Kugel am Gelenk
		stream << "sphere{<" << nodeIter->getX() << "," << nodeIter->getY()
				<< "," << nodeIter->getZ() << ">,r}" << endl;

		Neighbour * neighIter = nodeIter->getNeighbours();
		while (neighIter) {

			//  nur wenn nachbar und knoten nicht über die grenze verbunden sind...
			if (list->isPeriodic()
					&& nodeIter->euklidian(neighIter->getNode())
							< list->lengthX() / 2) {
				// Zylinder von A nach B TODO: den zurück nicht...
				stream << "cylinder{<" << nodeIter->getX() << ","
						<< nodeIter->getY() << "," << nodeIter->getZ() << ">,<"
						<< neighIter->getNode()->getX() << ","
						<< neighIter->getNode()->getY() << ","
						<< neighIter->getNode()->getZ() << ">,r}" << endl;

			} else if (!list->isPeriodic()) {
				// Zylinder von A nach B TODO: den zurück nicht...
				stream << "cylinder{<" << nodeIter->getX() << ","
						<< nodeIter->getY() << "," << nodeIter->getZ() << ">,<"
						<< neighIter->getNode()->getX() << ","
						<< neighIter->getNode()->getY() << ","
						<< neighIter->getNode()->getZ() << ">,r}" << endl;
			}

			neighIter = neighIter->getNext();
		}

		data.push_back(stream.str());

		nodeIter = nodeIter->getNext();
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
*/
/**
 * Schreibt die Koordinaten der Knoten mit einer Nummerierung davor, für voro++.
 *//*
void writePointsVoro(nodelist * list, const char outfileName[]) {

	vector<string> data;
	int ctr = 0;
	for (vector<node * >::iterator it = list->begin(); it != list->end(); ++it) {
	node * n = it;
		std::ostringstream strs;
		strs << ctr << "\t" << n->getPosition().toString("", "\t", "");

		data.push_back(strs.str());

		// weiter gehts
		ctr++;
	}
	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Daten schreiben
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i] << endl;
	}

	cout << "Muster in " << outfileName << " geschrieben." << endl;
}*/

/**
 * Test für voro++
 *//*
void testVoro() {
	// Set up constants for the container geometry
	const double x_min = -5, x_max = 5;
	const double y_min = -5, y_max = 5;
	const double z_min = -5, z_max = 5;

	// Set up the number of blocks that the container is divided into
	const int n_x = 6, n_y = 6, n_z = 6;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	voro::container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z,
			true, true, true, 8);

	//Randomly add particles into the container
	con.import("./data/voro++.dat");

	con.print_custom("%i %s", "./data/voro++_neighbours.dat");

	// Save the Voronoi network of all the particles to text files
	// in gnuplot and POV-Ray formats
	con.draw_cells_gnuplot("./data/voro++.gnu");
	con.draw_cells_pov("./data/voro++_c.pov");

	// Output the particles in POV-Ray format
	con.draw_particles_pov("./data/voro++_p.pov");
}

void compareLists(nodelist * listA, nodelist * listB) {
	cout << "Statistik Liste 1:" << endl << listA->listStats();
	cout << "Statistik Liste 2:" << endl << listB->listStats();

	// Mittelpunkte auf den von listA gleichsetzen
	coordinate shifter;
	listB->shiftList(
			(listB->getMins() + listB->getMaxs()
					- (listA->getMins() + listA->getMaxs())) / 2);

	// Dichte anpassen
	listB->setDensity(listA->getDensity());

	cout << "Statistik Liste 1:" << endl << listA->listStats();
	cout << "Statistik Liste 2:" << endl << listB->listStats();
}

void benchmark() {
	double pos1[] = { M_PI, 3.0, 5.0 };
	double pos2[] = { 5.0, 17.0, 1.0 };

	coordinate coord1(pos1, 2);
	coordinate coord2(pos1, 2);
	coordinate coord3(pos2, 2);

	cout << "1 " << coord1 << endl;
	cout << "2 " << coord2 << endl;
	cout << "3 " << coord3 << endl;

	clock_t t;
	t = clock();
	for (unsigned i = 0; i < 1000; i++) {
		(coord1 == coord2);
		(coord1 == coord3);
		(coord1 != coord2);
		(coord1 != coord3);
	}
}*/

/**
 * Hier wird ausgeführt was gewählt wurde.
 */
int main(int argc, char *argv[]) {
	cout << "Es gilt also ein Punktmuster zu charakterisieren. Also los!"
			<< endl;

	//vector<nodelist * > lists; // TODO
	// Muster einlesen
	nodelist * list = readfile(argv[1], argv[2]);
	nodelist * list2 = readfile(argv[3], argv[4]);

	int option = -1;
	while (option != 0) {
		// Optionen anzeigen und wählen lassen
		option = gui();

		switch (option) {
		case 0:
			cout << "Beende." << endl;
			break;
		case 1:
			list->neighbourDistribution();
			break;/*
		case 2:
			lengthDistribution(list,
					"./data/statistics/lenghtDistribution.dat");
			break;
		case 3:
			angleDistribution(list, "./data/statistics/angleDistribution.dat");
			break;
		case 4:
			hyperuniformity(list, "./data/statistics/hyperuniformity.dat");
			break;
		case 5:
			list->display();
			break;
		case 6:
			cout << "Liste 1: " << list->listStats();
			cout << "Liste 2: " << list2->listStats();
			break;
		case 7:
			gnuplotPattern(list);
			break;
		case 8:
			writePOV(list, "./data/staebe.pov");
			break;
		case 9:
			writePointsVoro(list, "./data/voro++.dat");
			break;
		case 10:
			testVoro();
			break;
		case 11:
			list2 = readfileGui();
			break;
		case 12:
			compareLists(list, list2);
			break;*/
		default:
			cout << "Das gibts leider (noch) nicht." << endl;
		}

		if (option != 0) {
			cout << "Und nun?" << endl;
		}
	}
	return 0;
}
