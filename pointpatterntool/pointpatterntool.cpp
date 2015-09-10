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
 * Reads two files in. The first contains coordinates in R^3, the second neighbours of those. Example:
 * nodeFile:	neighbourFile:
 * 1 2 3		2 3 4
 * 1 2 3		3 4 5
 * 2 3 4		1 2 3
 * 2 3 4		3 4 5
 * 3 4 5		1 2 3
 * 3 4 5		2 3 4
 */
nodelist* readfile(const char* nodes, const char* neighbours, bool periodic,
		string name) {
	cout << "Read nodesfile " << nodes << " and neighboursfile " << neighbours
			<< "." << endl;

	// input filestream of both files
	ifstream nodeFile, neighbourFile;
	nodeFile.open(nodes, std::ifstream::in);
	neighbourFile.open(neighbours, std::ifstream::in);

	// datastructure of the list
	nodelist* list = new nodelist(periodic, name);

	if (periodic) {
		cout << "Assume a box of dimensions (-5,5)^3..." << endl;
		list->setMins(coordinate(-5, -5, -5));
		list->setMaxs(coordinate(5, 5, 5));
	}

	// position of the node
	double x, y, z;
	// position of the neighbour
	double nx, ny, nz;

	// boundaries for nonperiodic patterns
	double inf = numeric_limits<double>::infinity();
	double minX = inf, minY = inf, minZ = inf, maxX = -inf, maxY = -inf, maxZ =
			-inf;

	// actual readig of the files
	while ((nodeFile >> x >> y >> z) && (neighbourFile >> nx >> ny >> nz)) {

		//**// FIXME in this loop, needed memory increases exorbitantly
		// If there is no node at the read position...
		node* n = list->getAt(coordinate(x, y, z));
		if (!n) {
			// ... add one at the position
			n = new node(list, x, y, z);
			list->push_back(n);
		}

		node* neighbour = list->getAt(coordinate(nx, ny, nz));
		// If there is no neighbour at the read position of the second file...
		if (!neighbour) {
			// ... add one at the position
			neighbour = new node(list, nx, ny, nz);
			list->push_back(neighbour);
		}
		//**//

		// set neighbourhood
		n->addNeighbour(neighbour);
		neighbour->addNeighbour(n);

		// evaluate boundaries for nonperiodic patterns
		if (!periodic) {
			minX = min(minX, x);
			minY = min(minY, y);
			minZ = min(minZ, z);

			maxX = max(maxX, x);
			maxY = max(maxY, y);
			maxZ = max(maxZ, z);
		}

	}
	// close files
	nodeFile.close();
	neighbourFile.close();

	// set boundaries for nonperiodic patterns and evaluate edgenodes
	if (!periodic) {
		list->setMins(coordinate(minX, minY, minZ));
		list->setMaxs(coordinate(maxX, maxY, maxZ));

		// adaptable parameter
		double characteristicLength = stats(list->lengthDistribution(false))[1]
				* 1.5;
		cout << "Nodes farther away of the edge than " << characteristicLength
				<< " are declared to be edgenodes." << endl;
		for (vector<node*>::iterator n = list->begin(); n != list->end(); ++n) {
			(*n)->setEdgenode(characteristicLength);
		}
	}

	/*
	 list->scaleList(1. / 50.);
	 list->shiftList(coordinate(-5,-5,-5));
	 */

	cout << "Pattern read and saved. Statistics:" << endl
			<< list->listStats() << endl;

	cout << "Needs " << list->capacity()*sizeof(node) + sizeof(list) << " bytes of memory." << endl;

	list =NULL;
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
					//if (!(*neighIter)->isEdgenode()) {
					//if (!(*nodeIter)->isEdgenode()) {
					data.push_back(*(*nodeIter)->getPosition());
					data.push_back(*(*neighIter)->getPosition());

					data.push_back(coordinate(nan(""), nan(""), nan("")));

					//}
					//}
				}
			}
		}
		datas.push_back(data);
	}

	// Plotten
	plot3D(datas);
}

/**
 * Vergleich zwei Listen.
 */
void compareLists(vector<nodelist*>& lists) {
	// Skalierung anpassen
	// TODO: evtl guiabfrage
	///////////////////////////////////////////////////////////////////////
	lists[1]->scaleList(1. / 50.);

	// Schwerpunkt grob anpassen
	lists[1]->shiftList(lists[1]->getMid() * -1);

	// Differenzvektoren bilden
	vector<coordinate> diffs;
	//coordinate diff;
	// Schwerpunkt der Differenzen
	coordinate centreOfMass(0, 0, 0);

	unsigned comparisons = min(lists[0]->size(), lists[1]->size());
	for (unsigned i = 0; i < comparisons; i++) {
		if (!(*lists[0])[i]->isEdgenode()) {

			coordinate* vec1 = (*lists[0])[i]->getPosition();

			// Vektor mit kleinstem Abstand in zweiter liste finden
			coordinate* vec2 = (*lists[1])[i]->getPosition();
			coordinate diff = *vec1 - *vec2;
			for (nodelist::iterator n = lists[1]->begin(); n != lists[1]->end();
					++n) {
				if ((*vec1 - *((*n)->getPosition())).lengthSqr()
						< diff.lengthSqr()) {
					vec2 = (*n)->getPosition();
					diff = *vec1 - *vec2;
				}
			}

			// Differenzvektor merken
			diffs.push_back(diff);
			// zum Schwerpunkt berechnen
			centreOfMass += diff;
		}
	}

	// Schwerpunkt in die Mitte setzen
	centreOfMass /= diffs.size();
	// ZWEITE Liste um POSITIVEN Schwerpunktsvektor verschieben
	lists[1]->shiftList(centreOfMass);

	// diffs anpassen, Längen berechnen
	vector<double> diffLengths;
	for (vector<coordinate>::iterator diffIter = diffs.begin();
			diffIter != diffs.end(); ++diffIter) {
		(*diffIter) -= centreOfMass;
		diffLengths.push_back((*diffIter).length());
	}

	// Diffdaten plotbar machen
	vector<vector<coordinate> > plotData;
	plotData.push_back(diffs);

	vector<vector<double> > histData;
	histData.push_back(diffLengths);

	vector<string> title;

	title.push_back("Differenz");
	/*
	 plotHist(histData, 0, 0.1, 50, title, "Differenzlängen");
	 // Punktewolke
	 plot3D(plotData, "x", "y", "z", "w p ls 7");
	 // Angepasste Muster
	 gnuplotPattern(lists);
	 */
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
	cout << "Pointpatterns are to be characterized. Here we go!" << endl;

	// Muster einlesen
	vector<nodelist*> lists;

	// keine Argumente übergeben -> generiere Zufallsmuster und Diamantmuster
	if (argc == 1) {
		lists.push_back(new nodelist(1, false));
		lists.push_back(new nodelist(2, false));
	}
	for (int i = 1; i < argc; i += 3) {
		string name = "Pattern";
		if (i == 1) {
			name = "HPU Chi 4C Chi 0.13 NP=1000 UC";
		} else if (i == 4) {
			name = "Validation recovered pattern";
		}
		lists.push_back(
				readfile(argv[i], argv[i + 1], convert(argv[i + 2], false),
						name));
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
		case 1: {
			cout << "Zähle die Nachbarn jedes Knotens..." << endl;
			vector<vector<double> > data;
			vector<string> names;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				vector<double> neighbours = (*list)->neighbourDistribution();
				data.push_back(neighbours);

				names.push_back((*list)->getName());

				cout << "Nachbarnstatistik Muster " << (*list)->getName() << ":"
						<< endl;
				cout << statsAsString(stats(neighbours));
			}

			plotHist(data, 0, 10, 10, names, "Nachbarn");

			break;
		}
		case 2: {
			cout << "Bestimme die Längen zwischen benachbarten Punkten..."
					<< endl;
			vector<vector<double> > data;
			vector<string> names;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				vector<double> lengths = (*list)->lengthDistribution();
				data.push_back(lengths);

				names.push_back((*list)->getName());

				cout << "Längenstatistik Muster " << (*list)->getName() << ":"
						<< endl;
				cout << statsAsString(stats(lengths));
			}
			plotHist(data, 0, 1, 50, names, "Längen");
			break;
		}
		case 3: {
			cout << "Bestimme die Winkel zwischen benachbarten Punkten..."
					<< endl;
			vector<vector<double> > data;
			vector<string> names;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				vector<double> angles = (*list)->angleDistribution();
				data.push_back(angles);

				names.push_back((*list)->getName());

				cout << "Winkelstatistik Muster " << (*list)->getName() << ":"
						<< endl;
				cout << statsAsString(stats(angles));
			}
			plotHist(data, 0, 180, 180, names, "Winkel");
			break;
		}
		case 4: {
			cout << "Bestimme den Grad der Hyperuniformity für " << endl;
			vector<vector<vector<double> > > variances;
			vector<string> names;
			double xMax = 0;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				cout << char(9) << (*list)->getName() << endl;
				vector<vector<double> > variance;
				(*list)->hyperuniformity(variance);
				variances.push_back(variance);
				names.push_back((*list)->getName());
				xMax = std::max(xMax, variance[variance.size() - 1][0]);
			}
			plotHyperuniformity(variances, xMax, names, "Radius R",
					"Variance  {/Symbol s}^2(R)");
			break;
		}
		case 5: {
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->display();
			}
			break;
		}
		case 6: {
			ctr = 0;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				cout << "Liste " << ctr << ": " << (*list)->getName() << endl
						<< char(9) << (*list)->listStats();
				ctr++;
			}
			break;
		}
		case 7: {
			gnuplotPattern(lists);
			break;
		}
		case 8: {
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->writePOV();
			}
			break;
		}

		case 9: {
			compareLists(lists);
			break;
		}
		case 10: {
			ctr = 0;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				double factor = (*list)->normalize();
				cout << "Muster " << ctr << "normalisiert (skaliert mit Faktor "
						<< factor << ").";
				ctr++;
			}
			break;
		}
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
