/*
 * pointpatterntool.cpp
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

#include <iostream>
#include <fstream>

#include "nodelist.hpp"

using namespace std;

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

	// Listenkopf
	NodeHead * list = new NodeHead();

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
	cout << "1 - Vert. der Anzahl der Nachbarn." << char(9) << "5 - " << endl;
	cout << "2 - Vert. der Abstände zu den Nachbarn." << char(9) << "6 - "
			<< endl;
	cout << "3 - Vert. der Winkel der Nachbarn." << char(9) << "7 - " << endl;
	cout << "4 - Grad der Hyperuniformity." << char(9) << char(9) << "8 - "
			<< endl << endl;
	cout << "9 - Liste darstellen" << char(9) << "0 - Nichts." << endl;
	cout << "Was möchtest du über die Punkte wissen? (1-9, default: 0) >> ";

	int option = 0;

	// Wenn das nächste Zeichen ein newline ist...
	if (cin.peek() == '\n') {
		// ...defaultwert setzen
		option = 0;
		// wenn kein int eingegeben wurde...
	} else if (!(cin >> option)) {
		// ...Fehler schmeißen
		cout << "Das war keine Zahl!" << endl;
	}

	return option;
}

/**
 * Schreibt die Anzahl der Nachbarn jedes Knotens untereinander in eine Datei und TODO führt ein gnuplot-script aus das die Daten als pdf plotet.
 */
void neighbourDistribution(NodeHead * list, const char outfileName[]) {
	cout << "Zähle die Nachbarn jedes Knotens..." << endl;

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Iterationsvariablen
	Node * nodeIter = list->getFirst();
	Neighbour * neighIter;
	int counter;

	// Iteration über jeden Knoten
	while (nodeIter) {
		counter = 0;
		neighIter = nodeIter->getNeighbours();
		// Nachbarn zählen
		while (neighIter) {
			counter++;
			neighIter = neighIter->getNextNeighbour();
		}
		// Wert schreiben
		outfile << counter << endl;
		nodeIter = nodeIter->getNext();
	}

	// TODO: gnuplotscript hier ausführen

	cout << "Nachbarsdaten in " << outfileName << " geschrieben." << endl;
}

/**
 * Schreibt alle Distanzen zwischen den Knoten und ihrer Nachbarn in eine Datei. TODO: gnuplotscript
 */
void lengthDistribution(NodeHead * list, const char outfileName[]) {
	cout << "Bestimme die Längen zwischen benachbarten Punkten..." << endl;

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Iterationsvariablen
	Node * nodeIter = list->getFirst();
	Neighbour * neighIter;

	// Iteration über jeden Knoten
	while (nodeIter) {
		neighIter = nodeIter->getNeighbours();
		// über Nachbarn iterieren
		while (neighIter) {
			// Wert schreiben
			outfile << nodeIter->euklidian(neighIter->getNode()) << endl;
			// weiter gehts
			neighIter = neighIter->getNextNeighbour();
		}
		nodeIter = nodeIter->getNext();
	}

	// TODO: periodische Randbedingungen beachten
	// TODO: doppelte Längen ausschließen (nach berechnung sortieren und jedes zweite Element schreiben?)
	// TODO: gnuplotscript hier ausführen

	cout << "Längendaten in " << outfileName << " geschrieben." << endl;
}

/**
 * Schreibt die Winkel zwischen den Nachbarn in eine Datei.
 */
void angleDistribution(NodeHead * list, const char outfileName[]) {
	cout << "Bestimme die Winkel zwischen benachbarten Punkten..." << endl;

	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Iterationsvariablen
	Node * nodeIter = list->getFirst();
	Neighbour * neighIter;

	// Iteration über jeden Knoten
	while (nodeIter) {
		neighIter = nodeIter->getNeighbours();
		// über Nachbarn iterieren
		while (neighIter) {
			while (neighIter->getNextNeighbour()) {
				// Wert schreiben
				outfile
						<< (180 / 3.141)
								* nodeIter->angle(neighIter->getNode(),
										neighIter->getNextNeighbour()->getNode())
						<< endl;
				neighIter = neighIter->getNextNeighbour();
			}

			// weiter gehts
			neighIter = neighIter->getNextNeighbour();
		}
		nodeIter = nodeIter->getNext();
	}

	// TODO: periodische Randbedingungen beachten
	// TODO: doppelte Winkel ausschließen (nach berechnung sortieren und jedes zweite Element schreiben?)
	// TODO: gnuplotscript hier ausführen

	cout << "Winkeldaten in " << outfileName << " geschrieben." << endl;
}

/**
 * Hier wird ausgeführt was gewählt wurde.
 */
int main(int argc, char *argv[]) {
	cout << "Es giilt also ein Punktmuster zu charakterisieren. Also los!"
			<< endl;
	NodeHead * list = readfile(argv[1], argv[2]);
	//list->display();
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
