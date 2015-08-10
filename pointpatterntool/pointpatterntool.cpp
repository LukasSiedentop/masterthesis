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

#include <voro++/voro++.hh>

#include <boost/tuple/tuple.hpp>
#include <boost/progress.hpp>

#include "templates.cpp"
#include "nodelist.hpp"

using namespace std;

/**
 * Plottet das Punktmuster.
 */
void plotPattern(NodeHead * list) {
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

	Gnuplot gp;
	gp << "reset\n";

	gp << "min = " << list->getMinX() - 0.5 << "\n";
	gp << "max = " << list->lengthX() + list->getMinX()+0.5 << "\n";

	gp << "set xlabel 'x'\n";
	gp << "set ylabel 'y'\n";
	gp << "set zlabel 'z'\n";

	gp << "set xrange [min:max]\n";
	gp << "set yrange [min:max]\n";
	gp << "set zrange [min:max]\n";

	// z-Achsen Offset ausschalten
	gp << "set ticslevel 0\n";
	gp << "set tics out nomirror\n";

	gp << "set xtics min,1,max\n";
	gp << "set ytics min,1,max\n";
	gp << "set ztics min,1,max\n";

	gp << "set view equal xyz\n";

	//gp << "set datafile missing 'nan'\n"; u ($1):($2):($3)
	gp << "splot '-' w l lc rgb'blue' notitle \n";
	gp.send1d(data);

	cout << "Weiter mit Enter." << endl;
	cin.get();
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
NodeHead * readfile(const char * nodes, const char * neighbours) {
	cout << "Lese Knotendatei " << nodes << " und Nachbardatei " << neighbours
			<< " ein..." << endl;

	// input filestream beider Dateien
	ifstream nodeFile, neighbourFile;
	nodeFile.open(nodes, std::ifstream::in);
	neighbourFile.open(neighbours, std::ifstream::in);

	cout << "Ist das Muster periodisch? (0,1; default: 0) >> ";

	bool periodic = input(0);

	// Listenkopf
	NodeHead * list = new NodeHead(periodic);

	if (periodic) {
		cout << "Nehme eine Box von (-5,5)^3 an..." << endl;
		list->setMins(-5, -5, -5);
		list->setMaxs(5, 5, 5);
	}

	// Position des Punktes
	double x, y, z;
	// Position des Nachbarn
	double nx, ny, nz;

	// einlesen der Dateien
	while ((nodeFile >> x >> y >> z) && (neighbourFile >> nx >> ny >> nz)) {

		Node * node = list->getAt(x, y, z);
		// Wenn an der Knoten Position noch kein Knoten existiert...
		if (!node) {
			// ... füge ihn der Liste hinzu
			node = new Node(x, y, z);
			list->add(node);
		}

		Node * neighbour = list->getAt(nx, ny, nz);
		// Wenn an der Nachbarknoten Position noch kein Knoten existiert...
		if (!neighbour) {
			// ... füge ihn der Liste hinzu
			neighbour = new Node(nx, ny, nz);
			list->add(neighbour);
		}

		// die Nachbarschaft setzen
		node->addNeighbour(neighbour);
		neighbour->addNeighbour(node);

	}

	cout << "Dateien eingelesen und Liste mit " << list->length()
			<< " Knoten erstellt." << endl;
	return list;
}

/**
 * Schreibt die Optionen auf die Koinsole und wartet auf einen Input.
 */
int gui() {
	cout << char(9) << "Statistiken" << char(9) << char(9) << char(9) << char(9)
			<< "Sortier- und Arrangieroptionen" << endl;
	cout << "1 - Vert. der Anzahl der Nachbarn." << char(9) << "10 - " << endl;
	cout << "2 - Vert. der Abstände zu den Nachbarn." << char(9) << "11 - "
			<< endl;
	cout << "3 - Vert. der Winkel der Nachbarn." << char(9) << "12 - " << endl;
	cout << "4 - Vert. der Voronoizellengröße." << char(9) << char(9) << "13 - "
			<< endl;
	cout << "5 - Grad der Hyperuniformity." << char(9) << char(9) << "14 - "
			<< endl;
	cout << "6 - Voro++ ausprobieren" << char(9) << "7 - Muster schreiben (voro++)" << char(9) << "8 - Muster plotten"
			<< char(9) << "9 - Liste darstellen" << char(9) << "0 - Nichts."
			<< endl;
	cout
			<< "Was möchtest du über die Punkte wissen? (0-14; default: 0) >> ";

	int option = 0;

	return input(option);
}

// TODO: Methoden die NodeHead * list brauchen in nodehead?!

/**
 * Schreibt die Anzahl der Nachbarn jedes Knotens untereinander in eine Datei.
 */
void neighbourDistribution(NodeHead * list, const char outfileName[]) {
	cout << "Zähle die Nachbarn jedes Knotens..." << endl;

	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Iterationsvariablen
	Node * nodeIter = list->getFirst();
	Neighbour * neighIter;
	int counter;

	// Iteration über jeden Knoten
	while (nodeIter) {
		counter = 0;
		neighIter = nodeIter->getNeighbours();
		// Nachbarn zählen
		// TODO: wenn nicht periodisch nur nachbarn die nicht am Rand sind.
		while (neighIter) {
			counter++;
			neighIter = neighIter->getNext();
		}
		// Wert schreiben
		data.push_back(counter);
		// weiter gehts
		nodeIter = nodeIter->getNext();
	}

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Daten schreiben
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i] << endl;
	}

	// Statistiken
	cout << "Nachbarstatistik:" << endl;
	stats(data);

	// Daten ploten
	plotHist(data, 0, 10, 10, "Anzahl Nächster Nachbarn");

	cout << "Nachbarsdaten in " << outfileName << " geschrieben." << endl;
}

/**
 * Schreibt alle Distanzen zwischen den Knoten und ihrer Nachbarn in eine Datei.
 */
void lengthDistribution(NodeHead * list, const char outfileName[]) {
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

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// ...und jeden zweiten schreiben.
	for (unsigned int i = 0; i < data.size(); i += 2) {
		outfile << data[i] << endl;
		halfData.push_back(data[i]);
	}

	// Statistiken
	cout << "Längenstatistik:" << endl;
	stats(halfData);

	// Daten ploten
	plotHist(halfData, 0.6, 1.05, 25, "Länge");

	cout << "Längendaten in " << outfileName << " geschrieben." << endl;
}

/**
 * Schreibt die Winkel zwischen den Nachbarn in eine Datei.
 */
void angleDistribution(NodeHead * list, const char outfileName[]) {
	cout << "Bestimme die Winkel zwischen benachbarten Punkten..." << endl;

	// Datenstruktur um die Werte aufzunehmen
	vector<double> data;

	// Fortschrittsbalken
	boost::progress_display show_progress(list->length());

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

		// Fortschritt
		++show_progress;
	}

	// Daten sortieren...
	data = mergeSort(data);
	vector<double> halfData;

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// ...und jeden zweiten schreiben.
	for (unsigned int i = 0; i < data.size(); i += 2) {
		halfData.push_back(data[i]);
		outfile << data[i] << endl;
	}

	// Statistiken
	cout << "Winkelstatistik:" << endl;
	stats(halfData);

	// Daten ploten
	plotHist(halfData, 0, 180, 180, "Winkel");

	cout << "Winkeldaten in " << outfileName << " geschrieben." << endl;
}

/**
 * Errechnet die Hyperuniformität des Musters.
 */
void hyperuniformity(NodeHead * list, const char outfileName[]) {

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
	outfile << "Radius" << char(9) << "Varianz" << endl;
	for (int i = 0; i < floor(rMax / dr); i++) {
		outfile << variance[i][0] << char(9) << variance[i][1] << endl;
	}

	// Varianz über Radius plotten
	plot2D(variance);

}

/**
 * Schreibt die Koordinaten der Knoten mit einer Nummerierung davor, für voro++.
 */
void writePoints(NodeHead * list, const char outfileName[]) {

	vector<string> data;
	int ctr = 0;
	Node * nodeIter = list->getFirst();
	while (nodeIter) {
		std::ostringstream strs;
		strs << ctr << "\t" << nodeIter->getX() << "\t" << nodeIter->getY()
				<< "\t" << nodeIter->getZ();

		data.push_back(strs.str());

		// weiter gehts
		ctr++;
		nodeIter = nodeIter->getNext();
	}

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Daten schreiben
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i] << endl;
	}

	cout << "Muster in " << outfileName << " geschrieben." << endl;
}

/**
 * Test für voro++
 */
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

/**
 * Hier wird ausgeführt was gewählt wurde.
 */
int main(int argc, char *argv[]) {

	cout << "Es gilt also ein Punktmuster zu charakterisieren. Also los!"
			<< endl;
	// Muster einlesen
	NodeHead * list = readfile(argv[1], argv[2]);

	int option = -1;
	while (option != 0) {
		// Optionen anzeigen und wählen lassen
		option = gui();

		switch (option) {
		case 0:
			cout << "Beende." << endl;
			break;
		case 1:
			neighbourDistribution(list,
					"./data/statistics/neighbourDistribution.dat");
			break;
		case 2:
			lengthDistribution(list,
					"./data/statistics/lenghtDistribution.dat");
			break;
		case 3:
			angleDistribution(list, "./data/statistics/angleDistribution.dat");
			break;
		case 5:
			hyperuniformity(list, "./data/statistics/hyperuniformity.dat");
			break;
		case 6:
			testVoro();
			break;
		case 7:
			writePoints(list, "./data/voro++.dat");
			break;
		case 8:
			plotPattern(list);
			break;
		case 9:
			list->display();
			break;
		default:
			cout << "Das gibts leider (noch) nicht." << endl;
		}

		if (option != 0) {
			cout << "Und nun?" << endl;
		}
	}

	return 0;
}
