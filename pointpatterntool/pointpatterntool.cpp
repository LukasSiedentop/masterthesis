/*
 * pointpatterntool.cpp
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

/*
 int main()
 {
 int* k; // Pointer to int
 int x =25; // int
 k= &x;
 cout<<k<<endl;     // Address of Operator
 cout<<*k<<endl;    // The Dereferencing Operator
 return 0;
 }





 func(int& a){
 do sth with a
 } -> die Referenz wird übergeben -> a wird selbst verändert


 new -> auf den heap
 nichts -> auf den stack

 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <ctime>

//#include <voro++/voro++.hh>

#include <boost/tuple/tuple.hpp>

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

	cout << "9 - Muster vergleichen" << endl;
	cout
			<< "10 - Muster normalisieren (Dichte = 1, Schwerpunkt in den Ursprung verschieben)"
			<< endl;

	cout << "0 - Nichts." << endl;
	cout << "---------------------------------------------------" << endl;
	cout << "Was möchtest du über die Punkte wissen? (0-10; default: 0) << ";

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
nodelist* readfile(const char* nodes, const char* neighbours, bool periodic) {
	cout << "Lese Knotendatei " << nodes << " und Nachbardatei " << neighbours
			<< " ein." << endl;

	// input filestream beider Dateien
	ifstream nodeFile, neighbourFile;
	nodeFile.open(nodes, std::ifstream::in);
	neighbourFile.open(neighbours, std::ifstream::in);

	// Listenkopf
	nodelist* list = new nodelist(periodic);

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
	double inf = numeric_limits<double>::infinity();
	double minX = inf, minY = inf, minZ = inf, maxX = -inf, maxY = -inf, maxZ =
			-inf;
	// einlesen der Dateien
	while ((nodeFile >> x >> y >> z) && (neighbourFile >> nx >> ny >> nz)) {
		node* n = list->getAt(coordinate(x, y, z));
		// Wenn an der Knoten Position noch kein Knoten existiert...
		if (!n) {
			// ... füge ihn der Liste hinzu
			n = new node(list, x, y, z);
			list->push_back(n);
		}

		node* neighbour = list->getAt(coordinate(nx, ny, nz));
		// Wenn an der Nachbarknoten Position noch kein Knoten existiert...
		if (!neighbour) {
			// ... füge ihn der Liste hinzu
			neighbour = new node(list, nx, ny, nz);
			list->push_back(neighbour);
		}

		// die Nachbarschaft setzen
		n->addNeighbour(neighbour);
		neighbour->addNeighbour(n);

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

	// Extremalwerte setzen und Randknoten bestimmen
	if (!periodic) {
		list->setMins(coordinate(minX, minY, minZ));
		list->setMaxs(coordinate(maxX, maxY, maxZ));

		// TODO: evtl doppelte oder halbe länge
		double characteristicLength = stats(list->lengthDistribution(false))[1]
				* 1.5;
		cout << "Knoten die vom Rand weiter weg sind als "
				<< characteristicLength << " werden als Randknoten deklariert."
				<< endl;
		for (vector<node*>::iterator n = list->begin(); n != list->end(); ++n) {
			(*n)->setEdgenode(characteristicLength);
		}
	}

	cout << (periodic ? "Periodisches " : "")
			<< "Muster eingelesen und Liste erstellt. Statistik:" << endl
			<< list->listStats() << endl;

	return list;
}

/**
 * Bereitet die Daten auf, das Gnuplot die Stäbe plotten kann.
 */
void gnuplotPattern(vector<nodelist*>& lists) {
	vector<vector<coordinate> > datas;

	// Listeniteration
	for (vector<nodelist*>::iterator list = lists.begin(); list != lists.end();
			++list) {
		vector<coordinate> data;
		// Knoteniteration
		for (vector<node*>::iterator nodeIter = (*list)->begin();
				nodeIter != (*list)->end(); ++nodeIter) {
			// Nachbarniteration
			for (vector<node*>::iterator neighIter =
					(*nodeIter)->getNeighbours()->begin();
					neighIter != (*nodeIter)->getNeighbours()->end();
					++neighIter) {

				//  nur wenn nachbar und knoten nicht über die grenze verbunden sind...
				if ((*list)->isPeriodic()
						&& (*nodeIter)->euklidian(*neighIter)
								< (*list)->getMaxFeatureSize()) {
					data.push_back(*(*nodeIter)->getPosition());
					data.push_back(*(*neighIter)->getPosition());

					data.push_back(coordinate(nan(""), nan(""), nan("")));

				} else if (!(*list)->isPeriodic()) {
					//if ((*neighIter)->isEdgenode()) {
					//if ((*nodeIter)->isEdgenode()) {
					data.push_back(*(*nodeIter)->getPosition());
					data.push_back(*(*neighIter)->getPosition());

					data.push_back(coordinate(nan(""), nan(""), nan("")));

				}
			}
		}
		datas.push_back(data);
	}

	// Plotten
	plot3D(datas);
}

void compareLists(vector<nodelist*>& lists) {
	// Skalierung anpassen
	// TODO: evtl guiabfrage
	lists[1]->scaleList(1. / 50.);

	// Schwerpunkt in die Mitte setzen
	lists[1]->shiftList(lists[1]->getMid() * -1);

	// Differenzvektoren bilden
	vector<vector<coordinate> > plotData;
	vector<coordinate> diffs;
	vector<double> diffLengths;
	coordinate diff;
	unsigned comparisons = min(lists[0]->size(), lists[1]->size());
	for (unsigned i = 0; i < comparisons; i++) {
		coordinate* vec1 = (*lists[0])[i]->getPosition();

		// Vektor mit kleinstem Abstand in zweiter liste finden
		coordinate* vec2 = (*lists[1])[i]->getPosition();
		diff = *vec1 - *vec2;

		for (nodelist::iterator n = lists[1]->begin(); n != lists[1]->end();
				++n) {
			if (!(*n)->isEdgenode()) {
				if ((*vec1 - *((*n)->getPosition())).lengthSqr()
						< diff.lengthSqr()) {
					vec2 = (*n)->getPosition();
					diff = *vec1 - *vec2;
				}
			}
		}

		diffs.push_back(diff);
		diffs.push_back(coordinate(0, 0, 0));
		diffs.push_back(coordinate(nan(""), nan(""), nan("")));

		diffLengths.push_back(diff.length());
	}

	plotData.push_back(diffs);

	// Histogramm der Längen der Differenzvektoren
	plotHist(diffLengths, 0, 0.1, 10, "Differenzlängen");
	plot3D(plotData, "x", "y", "z", "w dots");

	gnuplotPattern(lists);

}

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
 */
/*
 void benchmark() {
 double pos1[] = { M_PI, 3.0, 5.0 };
 double pos2[] = { 2, 2, 2 };
 double pos3[] = { 1.0, 1.0, 1.0 };

 coordinate coord1(pos1, 3);
 coordinate coord2(pos1, 3);
 coordinate coord3(pos2, 3);
 coordinate coord4(pos3, 3);

 cout << "1 " << coord1 << endl;
 cout << "2 " << coord2 << endl;
 cout << "3 " << coord3 << " lsqr: " << coord3.lengthSqr() << endl;
 cout << "4 " << coord4 << " lsqr: " << coord4.lengthSqr() << endl;
 cout << "3+4 " << coord4 + coord3 << " lsqr: "
 << (coord4 + coord3).lengthSqr() << endl;
 cout << "dist 3, 4 " << coord4.euklidian(coord3) << endl;
 cout << "3==4 " << (coord3 == coord4) << endl;
 cout << "4==3 " << (coord4 == coord3) << endl;
 cout << "3==3 " << (coord3 == coord3) << endl;

 clock_t t;
 t = clock();
 for (unsigned i = 0; i < 1000; i++) {
 (coord3 == coord4);
 }

 t = clock() - t;
 cout << "Hat " << t << " Clicks gedauert ("
 << (((float) t) / CLOCKS_PER_SEC) << "s)" << endl;
 }
 */

/**
 * Hier wird ausgeführt was gewählt wurde.
 */
int main(int argc, char* argv[]) {
	cout << "Es gilt Punktmuster zu charakterisieren. Also los!" << endl;

	// Muster einlesen
	vector<nodelist*> lists;
	for (int i = 1; i < argc; i += 3) {
		lists.push_back(
				readfile(argv[i], argv[i + 1], convert(argv[i + 2], false)));
	}

	int option = -1;
	while (option != 0) {
		//benchmark();
		int ctr;

		// Optionen anzeigen und wählen lassen
		option = gui();

		switch (option) {
		case 0:
			cout << "Beende." << endl;
			break;
		case 1:
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->neighbourDistribution();
			}
			break;
		case 2:
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->lengthDistribution();
			}
			break;
		case 3:
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->angleDistribution();
			}
			break;
		case 4:
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->hyperuniformity();
			}
			break;
		case 5:
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->display();
			}
			break;
		case 6:
			ctr = 0;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				cout << "Liste " << ctr << ": " << (*list)->listStats();
				ctr++;
			}
			break;
		case 7:
			gnuplotPattern(lists);
			break;
		case 8:
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->writePOV();
			}
			break;

		case 9:
			compareLists(lists);
			break;
		case 10:
			ctr = 0;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				double factor = (*list)->normalize();
				cout << "Muster " << ctr << "normalisiert (skaliert mit Faktor "
						<< factor << ").";
				ctr++;
			}
			break;
		default:
			cout << "Das gibts leider (noch) nicht." << endl;
		}

		if (option != 0) {
			cout << "Und nun?" << endl;
		}
	}

	// TODO: alles löschen
	for (vector<nodelist*>::iterator list = lists.begin(); list != lists.end();
			++list) {
		(*list)->~nodelist();
	}

	return 0;
}
