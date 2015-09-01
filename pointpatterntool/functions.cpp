/*
 * functions.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "functions.hpp"

void plot3D(vector<vector<coordinate> > datas, const char xlabel[],
		const char ylabel[], const char zlabel[]) {
	Gnuplot gp;
	gp << "reset\n";

	gp << "set xlabel '" << xlabel << "'\n";
	gp << "set ylabel '" << ylabel << "'\n";
	gp << "set zlabel '" << zlabel << "'\n";

	// z-Achsen Offset ausschalten
	gp << "set ticslevel 0\n";
	gp << "set tics out nomirror\n";

	gp << "set autoscale\n";

	gp << "set view equal xyz\n";

	// Farben
	vector<string> colors;
	colors.push_back("#1B9E77");
	colors.push_back("#D95F02");
	colors.push_back("#7570B3");
	colors.push_back("#E7298A");
	colors.push_back("#66A61E");
	colors.push_back("#E6AB02");
	colors.push_back("#A6761D");
	colors.push_back("#666666");

	// Plotstring bauen und an Gnuplot schicken
	stringstream plotstring;
	plotstring << "splot ";
	for (unsigned i = 0; i < datas.size(); i++) {
		plotstring << "'-' w l lt rgb '" << colors[i] << "' t 'Muster " << i << "'";
		if (i != datas.size() - 1) {
			plotstring << ",";
		}
	}
	gp << plotstring.str() << "\n";

	// Daten senden
	for (vector<vector<coordinate> >::iterator data = datas.begin();
			data != datas.end(); ++data) {

		gp.send1d(*data);
	}

	cout << "Weiter mit Enter." << endl;
	cin.get();
}
