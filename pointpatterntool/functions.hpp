/*
 * templates.hpp
 *
 *  Created on: 03.08.2015
 *      Author: Lukas Siedentop
 *
 * A Collection of useful methods like I/O, sorting etc that can be used frequently by methods from main.
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_

#include <vector>

#include <string>
#include <sstream>

#include <typeinfo>

#include <map>

#include "gnuplot-iostream.h"
#include "nodelist.hpp"
#include "coordinate.hpp"

using namespace std;

/* Statistiken */

/**
 * Berechnet die Statistiken eine Zahlenliste unter der Annahme, die Liste ist von klein nach groß Sortiert. TODO: Rückgabe evtl als map
 * Vektorposition Bedeutung
 * 0		Modus (TODO) https://de.wikipedia.org/wiki/Modus_(Statistik)
 * 1		Erwartungswert
 * 2		Varianz
 * 3		Schiefe
 * 4		Wölbung
 * 10		Anzahl
 * 11		Minimum
 * 12		Maximum
 * 13		Median
 * 14		variance Median
 */
vector<double> stats(vector<double> data);

/**
 * Gibt eine Liste mit Statistiken als menschenlesbaren String zurück.
 */
string statsAsString(const vector<double>& data, const char commentDelimeter[] =
		"\t");

/* I/O */

/*
 * Holt den Input Typgerecht. Von http://www.cplusplus.com/forum/articles/6046/
 * in muss mit dem defaultwert initialisiert sein, bei keiner Angabe wird dieser zurückgegeben.
 */
template<typename T>
T input(T in) {
	string input = "";

	while (true) {
		// füllt den input-string mit der Eingabe
		getline(cin, input);

		// setzt den Defaultwert bei einem leeren String (getline löscht das letzte \n!)
		if (input == "") {
			return in;
		}

		// Konvertieren von string nach T
		stringstream inputStream(input);
		if (inputStream >> in) {
			break;
		}
		cout << "Falscher Input! Nochmal." << endl;
	}

	return in;
}

/*
 * Konvertiert den gegebenen String in den Typ des gegebenen Arguments. Bei Fehler gibt kommt das Argument zurück.
 */
template<typename T>
T convert(string input, T type) {
	// Konvertieren von string nach T
	stringstream inputStream(input);
	if (inputStream >> type) {
		return type;
	}

	// Error
	cout << "Konnte string '" << input << "' nicht in " << typeid(type).name()
			<< " umwandeln, nehme " << type << "." << endl;
	return type;
}

/**
 * Schreibt das Histogramm in eine Datei mit gegebenem Namen, optional mit Statistiken und/oder Header
 */
template<typename T>
void writeHist(vector<T> data, bool includeStats, const char header[],
		const char outfileName[]) {
	// Outfile
	ofstream outfile;
	outfile.open(outfileName);

	// Statistiken
	if (includeStats) {
		outfile << statsAsString(stats(data), "# ") << endl;
	}

	// Kopfzeile
	outfile << header << endl;

	// Daten schreiben
	for (unsigned int i = 0; i < data.size(); i++) {
		outfile << data[i] << endl;
	}

	cout << "Histogramm in " << outfileName << " geschrieben." << endl;
}

// TODO: write 2D, 3D,

/* Sortiermethoden */

/**
 * Templates um einen Vektor zu (merge-) sortieren. Geändert von http://www.bogotobogo.com/Algorithms/mergesort.php, code A
 */
/**
 *
 */
template<typename T>
vector<T> merge(vector<T> left, vector<T> right) {
	vector<T> result;
	while ((int) left.size() > 0 || (int) right.size() > 0) {
		if ((int) left.size() > 0 && (int) right.size() > 0) {

			// Hier passiert der Vergleich
			if ((T) left.front() <= (T) right.front()) {
				result.push_back((T) left.front());
				left.erase(left.begin());
			} else {
				result.push_back((T) right.front());
				right.erase(right.begin());
			}

		} else if ((int) left.size() > 0) {
			for (int i = 0; i < (int) left.size(); i++)
				result.push_back(left[i]);
			break;
		} else if ((int) right.size() > 0) {
			for (int i = 0; i < (int) right.size(); i++)
				result.push_back(right[i]);
			break;
		}
	}
	return result;
}

template<typename T>
vector<T> mergeSort(vector<T> m) {
	if (m.size() <= 1)
		return m;

	vector<T> left, right, result;
	int middle = ((int) m.size() + 1) / 2;

	for (int i = 0; i < middle; i++) {
		left.push_back(m[i]);
	}

	for (int i = middle; i < (int) m.size(); i++) {
		right.push_back(m[i]);
	}

	left = mergeSort(left);
	right = mergeSort(right);
	result = merge(left, right);

	return result;
}

/* Plotmethoden Gnuplot */
// TODO: anzeigen UND svg speichern
/**
 * Plottet ein Histogramm der gegebenen (1D-)Daten mit der Binsize (max-min)/n und der x-Achsen-Beschriftung xlabel.
 */
template<typename T>
void plotHist(vector<T> data, double min, double max, int n,
		const char xlabel[] = "x") {
	Gnuplot gp;
	// Don't forget to put "\n" at the end of each line!

	gp << "reset\n";

	gp << "n=" << n << "\n";
	gp << "max=" << max << "\n";
	gp << "min=" << min << "\n";
	gp << "width=" << (max - min) / n << "\n";

	gp << "hist(x,width)=width*floor(x/width) #+width/2.0\n";

	gp << "set xrange [min:max]\n";
	gp << "set yrange [0:]\n";
	//gp << "set offset graph 0.05,0.05,0.05,0.0\n";
	// 10 x-tics
	gp << "set xtics min," << (max - min) / 10 << ",max\n";
	gp << "set boxwidth width*0.9\n";
	gp << "set style fill solid 0.5\n";
	gp << "set tics out nomirror\n";
	gp << "set xlabel '" << xlabel << "'\n";
	gp << "set ylabel 'Häufigkeit'\n";

	// '-' means read from stdin.  The send1d() function sends data to gnuplot's stdin.
	gp
			<< "plot '-' u (hist($1, width)):(1.0) w boxes smooth freq lc rgb'blue' notitle\n";
	gp.send1d(data);

	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	cout << "Weiter mit Enter." << endl;
	cin.get();
}

/**
 * Plottet die Varianz über den Radius. TODO: klären welches vector verwendet wird. ist das von boost besser? müsste immo std sein.
 */
template<typename T>
void plot2D(vector<vector<T> > data, double xMin, double dx, double xMax,
		double fitMax, const char xlabel[] = "x", const char ylabel[] = "y") {

	Gnuplot gp;
	gp << "reset\n";

	gp << "A = 1\n";
	gp << "f(x) = A*x**2\n";
	//gp << "f(x) = A*log(B*x)\n";
	gp << "fit [0:" << fitMax << "] f(x) '-' u 1:2 via A\n";
	gp.send1d(data);

	gp << "set key top left\n";
	gp << "set xrange [" << xMin << ":" << xMax << "]\n";
	gp << "set yrange [0:]\n";
	//gp << "set offset graph 0.05,0.05,0.05,0.0\n";
	// 10 x-tics
	gp << "set xtics " << xMin << ",1," << xMax << "\n";
	gp << "set tics out nomirror\n";
	gp << "set xlabel '" << xlabel << "'\n";
	gp << "set ylabel '" << ylabel << "'\n";

	// '-' means read from stdin.  The send1d() function sends data to gnuplot's stdin.
	gp
			<< "plot '-' u 1:2 ls 7 lc rgb'blue' t 'Numerische Daten', f(x) lc rgb'blue' t sprintf('Fit: f(R) = %.3fR^2',A)\n";
	gp.send1d(data);

	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Weiter mit Enter." << std::endl;
	std::cin.get();
}
;

/**
 * Plottet die gegebenen Daten im Format Spalten(Zeilen[3]) mit den gegebenen Grenzen.
 */
void plot3D(vector<vector<coordinate> > data, const char xlabel[] = "x",
		const char ylabel[] = "y", const char zlabel[] = "z",
		const char style[] = "w l");

#endif /* NODELIST_HPP_ */
