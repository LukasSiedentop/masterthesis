/*
 * functions.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "functions.hpp"


void plot3D(vector<coordinate > data, coordinate mins, coordinate maxs,
		const char xlabel[], const char ylabel[],
		const char zlabel[]) {
	Gnuplot gp;
	gp << "reset\n";

	gp << "set xlabel '" << xlabel << "'\n";
	gp << "set ylabel '" << ylabel << "'\n";
	gp << "set zlabel '" << zlabel << "'\n";

	gp << "set xrange [" << mins[0] - 0.5 << ":" << maxs[0] + 0.5 << "]\n";
	gp << "set yrange [" << mins[1] - 0.5 << ":" << maxs[1] + 0.5 << "]\n";
	gp << "set zrange [" << mins[2] - 0.5 << ":" << maxs[2] + 0.5 << "]\n";

	// z-Achsen Offset ausschalten
	gp << "set ticslevel 0\n";
	gp << "set tics out nomirror\n";

	gp << "set xtics " << mins[0] - 0.5 << ",1," << maxs[0] + 0.5 << "\n";
	gp << "set ytics " << mins[1] - 0.5 << ",1," << maxs[1] + 0.5 << "\n";
	gp << "set ztics " << mins[2] - 0.5 << ",1," << maxs[2] + 0.5 << "\n";

	gp << "set view equal xyz\n";

	//gp << "set datafile missing 'nan'\n"; u ($1):($2):($3)
	gp << "splot '-' w l lc rgb'blue' notitle \n";
	gp.send1d(data);

	cout << "Weiter mit Enter." << endl;
	cin.get();
}
