/*
 * pointpatterntool.cpp
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

/*
 * Comments should say "why", not "what".
 *
 *
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
	cout << "1 - number of neighbours histogram" << endl;
	cout << "2 - distance to neighbours histogram" << endl;
	cout << "3 - angle between neighbours histogram" << endl;
	cout << "4 - degree of hyperuniformity" << endl;

	cout << "5 - output pattern to console" << endl;
	cout << "6 - pattern statistics" << endl;
	cout << "7 - gnuplot pattern" << endl;
	cout << "8 - write POV-Ray file" << endl;

	cout << "9 - compare two patterns" << endl;
	cout
			<< "10 - normalize pattern (density of points= 1, shift center of mass to origin)"
			<< endl;

	cout << "11 - write MEEP dielectric" << endl;
	cout << "12 - write GWL file" << endl;
	cout << "13 - write MPB dielectric" << endl;
	cout << "14 - write ASCII file of extended pattern" << endl;

	cout << "15 - add random pattern" << endl;
	cout << "16 - add diamond pattern" << endl;

	cout << "17 - read glass file" << endl;

	cout << "0 - Nothing." << endl;
	cout << "---------------------------------------------------" << endl;
	cout << "What do you want to do? (0-17; default: 0) << ";

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
		std::string name) {
	std::cout << "Read nodesfile " << nodes << " and neighboursfile " << neighbours
			<< "." << std::endl;

	// input filestream of both files
	ifstream nodeFile, neighbourFile;
	nodeFile.open(nodes, std::ifstream::in);
	neighbourFile.open(neighbours, std::ifstream::in);

	// datastructure of the list
	nodelist* list = new nodelist(periodic, false, name);

	if (periodic) {
		std::cout << "Assume a box of dimensions (-5,5)^3..." << std::endl;
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

	// actual reading of the files
	while ((nodeFile >> x >> y >> z) && (neighbourFile >> nx >> ny >> nz)) {
		node* n = list->add(x, y, z);
		node* neighbour = list->add(nx, ny, nz);

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
		double characteristicLength = stats(list->lengthDistribution())[1]
				* 2;
		std::cout << "Nodes farther away of the edge than " << characteristicLength
				<< " are declared to be edgenodes." << std::endl;
		list->setEdgenodes(characteristicLength);
	}

	/*
	 list->scaleList(1. / 50.);
	 list->shiftList(coordinate(-5,-5,-5));
	 */

	std::cout << "Pattern read. Statistics:" << std::endl << list->listStats();

	//cout << "Needs " << list->capacity() * sizeof(node) + sizeof(list) << " bytes of memory." << endl;

	return list;
}

/**
 * Reads glass files as provided by Antonio Puertas, see /home/lukas/data/pointpatterns/glass/puertas/README.md
 */
nodelist* readglassfile(const char* nodes, std::string name) {
	std::cout << "Read nodesfile " << nodes << "." << std::endl;

	// input filestream of both files
	ifstream nodeFile;
	nodeFile.open(nodes, std::ifstream::in);

	// datastructure of the list
	nodelist* list = new nodelist(true, false, name);

	// calculation time (unimportant parameter) and boxlenght L
	double L, time;
	if (!(nodeFile >> time >> L)) {
		std::cout << "Error: head line is not in the format '  <time>  <sidelength>'!" << std::endl;
		return list;
	}

	list->setMins(coordinate((double) -L/2, (unsigned int) 3));
	list->setMaxs(coordinate((double) L/2, (unsigned int) 3));

	// position, velocity and radius of the particle
	double x, y, z, vx, vy, vz, r;

	// actual reading of the files
	while ((nodeFile >> x >> y >> z >> vx >> vy >> vz >> r)) {
		list->add(x, y, z);
	}

	// close files
	nodeFile.close();

	// normalize to 10x10x10 box
	std::cout << "Normalize from " << L << "^3 cube to 10^3 cube." << std::endl;
	list->scaleList(10/L);

	// define neighbours
	list->setNeighbours(4);


	std::cout << "Glass file read and neighbours constructed. Statistics:" << std::endl << list->listStats();

	//cout << "Needs " << list->capacity() * sizeof(node) + sizeof(list) << " bytes of memory." << endl;

	return list;
}

/**
 * Prepare data so that gnuplot can display the pattern.
 */
void gnuplotPattern(vector<nodelist*>& lists) {
	vector<vector<vector<double> > > datas;

	// listiteration
	for (vector<nodelist*>::iterator list = lists.begin(); list != lists.end();
			++list) {
		datas.push_back((*list)->getGnuplotMatrix());
	}

	// get names of lists
	std::vector<std::string> names;
	for (vector<nodelist*>::iterator list = lists.begin(); list != lists.end();
			++list) {
		names.push_back((*list)->getName());
	}

	plot3D(datas, names);
}

/**
 * Compares two lists. TODO: compare more than two with the first...
 */
void compareLists(vector<nodelist*>& lists) {

	// get names of lists
	vector<string> names;
	for (vector<nodelist*>::iterator list = lists.begin(); list != lists.end();
			++list) {
		names.push_back((*list)->getName());
	}

	// neighbourstats
	cout << "Count the neighbours of every node..." << endl;
	vector<vector<double> > neighData;
	for (vector<nodelist*>::iterator list = lists.begin();
			list != lists.end(); ++list) {
		vector<double> neighbours = (*list)->neighbourDistribution();
		neighData.push_back(neighbours);

		cout << "Neighbourstatistics of pattern '" << (*list)->getName()
				<< "':" << endl;
		cout << statsAsString(stats(neighbours));
	}

	plotHist(neighData, 0, 10, 10, names, "number of neighbours",
						"data/resolution-validation/neighbours_histogram");


	// length stats
	cout << "Evaluate the distance between neighbouring nodes..."
					<< endl;
	vector<vector<double> > lenData;
	for (vector<nodelist*>::iterator list = lists.begin();
			list != lists.end(); ++list) {
		vector<double> lengths = (*list)->lengthDistribution();
		lenData.push_back(lengths);

		cout << "Lengthstatistics of pattern " << (*list)->getName()
				<< ":" << endl;
		cout << statsAsString(stats(lengths));
	}
	plotHist(lenData, 0, 1, 50, names, "distance of neighbouring nodes",
					"data/resolution-validation/lengths_histogram");

	// angle stats
	cout << "Evaluate the angle between neighbouring nodes..." << endl;
	vector<vector<double> > angData;
	for (vector<nodelist*>::iterator list = lists.begin();
			list != lists.end(); ++list) {
		vector<double> angles = (*list)->angleDistribution();
		angData.push_back(angles);

		names.push_back((*list)->getName());
		cout << "Anglestatistics of pattern " << (*list)->getName()
				<< ":" << endl;
		cout << statsAsString(stats(angles));
	}
	plotHist(angData, 0, 180, 180, names,
			"angle between neighbouring nodes",
			"data/resolution-validation/angles_histogram");

	// TODO: evtl guiabfrage
	///////////////////////////////////////////////////////////////////////
	// adjust scaling and roughly the midpoint
	//lists[1]->scaleList(1. / 50.);
	//lists[1]->shiftList(lists[1]->getMid() * -1);

	vector<coordinate> diffs;
	// center of dass of difference vectors
	coordinate centreOfMass(0, 0, 0);

	unsigned comparisons = min(lists[0]->size(), lists[1]->size());
	for (unsigned i = 0; i < comparisons; i++) {

		if (!(*lists[0])[i].isEdgenode()) {

			coordinate vec1 = (*lists[0])[i].getPosition();

			// get closest node
			coordinate vec2 = (*lists[1])[i].getPosition();
			coordinate diff = vec1 - vec2;
			for (vector<node*>::iterator n = lists[1]->begin();
					n != lists[1]->end(); ++n) {
				if ((vec1 - ((*n)->getPosition())).lengthSqr()
						< diff.lengthSqr()) {
					vec2 = (*n)->getPosition();
					diff = vec1 - vec2;
				}
			}

			diffs.push_back(diff);
			centreOfMass += diff;
		}
	}

	// shift the SECOND list about the POSITIVE center of mass
	centreOfMass /= diffs.size();
	lists[1]->shiftList(centreOfMass);

	vector<vector<double> > plotDiffs;

	// shift diffs, calculate lengths
	vector<double> diffLengths;
	for (vector<coordinate>::iterator diffIter = diffs.begin();
			diffIter != diffs.end(); ++diffIter) {
		(*diffIter) -= centreOfMass;
		plotDiffs.push_back(*(*diffIter).getVector());
		diffLengths.push_back((*diffIter).length());
	}

	// make it plottable
	vector<vector<vector<double> > > plotData;
	plotData.push_back(plotDiffs);

	vector<vector<double> > histData;
	histData.push_back(diffLengths);

	vector<string> title;

	title.push_back(names[1]);

	plotHist(histData, 0, 0.1, 50, title, "difference vector lengths", "data/resolution-validation/difference-length");
	// Punktewolke
	plot3D(plotData, names, "x", "y", "z", "w p ls 7");
	// Angepasste Muster
	//gnuplotPattern(lists);

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
 coordinate coord1(M_PI, 3.0, 5.0);
 coordinate coord2(M_PI, 3.0, 5.0);
 coordinate coord3(2, 2, 2);
 coordinate coord4(1.0, 1.0, 1.0);

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
 }*/

/**
 * Hier wird ausgeführt was gewählt wurde.
 */
int main(int argc, char* argv[]) {
	cout << "Pointpatterns are to be characterized. Here we go!" << endl;

	vector<nodelist*> lists;

	// poisson
	lists.push_back(new nodelist(1, true));

	// No arguments given -> generate diamond and random point pattern
	/*
	if (argc == 1) {
		// poisson
		lists.push_back(new nodelist(1, true));
		// diamond
		lists.push_back(new nodelist(2, true));
	}
	*/



	// argument 1: nodes file, 2: neighbours files, 3: isperiodic, 4: name
	for (int i = 1; i < argc; i += 4) {

		lists.push_back(
				readfile(argv[i], argv[i + 1], convert(argv[i + 2], false),
						argv[i+3]));

/*
		// for resolution validation
		if (i==5) {
			// res600
			lists.back()->scaleList(1./50.);
			lists.back()->shiftList(coordinate(-5,-5,-5));
		} else if (i==9) {
			// res300
			lists.back()->scaleList(1./25.);
			lists.back()->shiftList(coordinate(-5,-5,-5));
		} else if (i==13) {
			// res150
			lists.back()->scaleList(1./12.5);
			lists.back()->shiftList(coordinate(-5,-5,-5));
		} else if (i==17) {
			// res120
			lists.back()->scaleList(1./10.);
			lists.back()->shiftList(coordinate(-5,-5,-5));
		} else if (i==21) {
			// res100
			lists.back()->scaleList(1./(100./12.));
			lists.back()->shiftList(coordinate(-5,-5,-5));
		} else if (i==25) {
			// res75
			lists.back()->scaleList(1./6.25);
			lists.back()->shiftList(coordinate(-5,-5,-5));
		}*/
	}

	// scale back to 10^3 cube
	//lists.back()->scaleList(1./(6.25));
	// shift origin to mid of cube
	//lists.back()->shiftList(coordinate(-5,-5,-5));

	// normalize agapornis
	//lists.back()->normalize();

	// functions TODO: https://de.wikipedia.org/wiki/Radiale_Verteilungsfunktion, Structure factor, beautify periodic boundary conditions, generation tool
	// write pattern to file, to disable double entries, making reading better...
	int option = -1;
	while (option != 0) {
		//benchmark();
		int ctr;

		// Optionen anzeigen und wählen lassen
		option = gui();

		switch (option) {
		case 0:
			cout << "Exit." << endl;
			break;
		case 1: {
			cout << "Count the neighbours of every node..." << endl;
			vector<vector<double> > data;
			vector<string> names;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				vector<double> neighbours = (*list)->neighbourDistribution();
				data.push_back(neighbours);

				names.push_back((*list)->getName());

				cout << "Neighbourstatistics of pattern '" << (*list)->getName()
						<< "':" << endl;
				cout << statsAsString(stats(neighbours));
			}

			plotHist(data, 0, 10, 10, names, "number of neighbours",
					"data/statistics/neighbours_histogram");

			break;
		}
		case 2: {
			cout << "Evaluate the distance between neighbouring nodes..."
					<< endl;
			vector<vector<double> > data;
			vector<string> names;

			double maxLength = 0;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				vector<double> lengths = (*list)->lengthDistribution();
				data.push_back(lengths);

				names.push_back((*list)->getName());

				cout << "Lengthstatistics of pattern " << (*list)->getName()
						<< ":" << endl;
				cout << statsAsString(stats(lengths));
				maxLength = max(maxLength,stats(lengths)[12]);
			}

			cout << "Maximum length: " << maxLength << endl;

			plotHist(data, 0, ceil(maxLength), 50, names, "distance of neighbouring nodes",
					"data/statistics/lengths_histogram");
			break;
		}
		case 3: {
			cout << "Evaluate the angle between neighbouring nodes..." << endl;
			vector<vector<double> > data;
			vector<string> names;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				vector<double> angles = (*list)->angleDistribution();
				data.push_back(angles);

				names.push_back((*list)->getName());

				cout << "Anglestatistics of pattern " << (*list)->getName()
						<< ":" << endl;
				cout << statsAsString(stats(angles));
			}
			plotHist(data, 0, 180, 180, names,
					"angle between neighbouring nodes",
					"data/statistics/angles_histogram");
			break;
		}
		case 4: {
			unsigned int nr = 50, n = 100;
			cout << "Number of spheres (default: 100) << ";
			input(n);
			cout << "Number of radii (default: 50) << ";
			input(nr);

			cout << "Determine the degree of hyperuniformity for " << endl;
			vector<vector<vector<double> > > variances;
			vector<std::string> names;
			double xMax = 0;
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				cout << char(9) << (*list)->getName() << endl;
				vector<vector<double> > variance;
				variance = (*list)->hyperuniformity(nr, n);
				variances.push_back(variance);
				names.push_back((*list)->getName());
				xMax = std::max(xMax, variance[variance.size() - 1][0]);
			}

			stringstream filename;
			filename << "data/statistics/hyperuniformity_n=" << n << "_nr="
					<< nr;

			plotHyperuniformity(variances, xMax, names, "radius R",
					"variance  {/Symbol s}^2(R)", filename.str());

			//writeHyperuniformity(variances, names, nr, n);

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
				cout << "List " << ctr << ": " << (*list)->getName() << endl
						<< (*list)->listStats();
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
				cout << "Pattern " << ctr << "normalised (scaled with factor "
						<< factor << ")." << endl;
				ctr++;
			}
			break;
		}

		case 11: {
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->writeMEEP();
			}
			break;
		}

		case 12: {
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->writeGWL();
			}
			break;
		}

		case 13: {
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->writeMPB();
			}
			break;
		}

		case 14: {
			for (vector<nodelist*>::iterator list = lists.begin();
					list != lists.end(); ++list) {
				(*list)->writeCoordinates();
			}
			break;
		}

		case 15: {
			lists.push_back(new nodelist(1, true));
			break;
		}

		case 16: {
			lists.push_back(new nodelist(2, true));
			break;
		}


		case 17: {
			cout << "Read glass file as provided by Antonio Puertas." << endl;
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic050/pos-eq-050-01.dat", "conf1Phi_c=0.50"));
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic050/pos-eq-050-02.dat", "conf2Phi_c=0.50"));
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic050/pos-eq-050-03.dat", "conf3Phi_c=0.50"));
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic050/pos-eq-050-04.dat", "conf4Phi_c=0.50"));
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic050/pos-eq-050-05.dat", "conf5Phi_c=0.50"));

			/*
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic053/pos-eq-053-01.dat", "Phi_c=0.53"));
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic055/pos-eq-055-01.dat", "Phi_c=0.55"));
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic057/pos-eq-057-01.dat", "Phi_c=0.57"));
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic058/pos-eq-058-01.dat", "Phi_c=0.58"));
			lists.push_back(readglassfile("/home/lukas/data/pointpatterns/raw/glass/puertas/phic0585/pos-eq-0585-01.dat", "Phi_c=0.585"));
			*/
			break;
		}

		default:
			cout << "This option does not (yet) exist, unfortunately." << endl;
		}

		if (option != 0) {
			cout << "What's next?" << endl;
		}
	}

	for (vector<nodelist*>::iterator list = lists.begin(); list != lists.end();
			++list) {
		(*list)->deleteEntries();
		delete (*list);
	}

	return 0;
}
