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
#include <algorithm>
#include <typeinfo>
#include <iostream>
#include <map>

#include "gnuplot-iostream.h"
#include "nodelist.hpp"
#include "coordinate.hpp"


//#include <CGAL/IO/Geomview_stream.h>

/* Statistiken */

/**
 * Calculates the statistics of a list of numbers with the asumption the list is sorted from small to big.
 * TODO: return map
 * position meaning
 * 0		mode (TODO) https://de.wikipedia.org/wiki/Modus_(Statistik)
 * 1		expectancy
 * 2		variance
 * 3		skewness
 * 4		kurtosis
 * 10		count
 * 11		minimum
 * 12		maximum
 * 13		median
 * 14		variance median
 */
std::vector<double> stats(std::vector<double> data);

/**
 * Gibt eine Liste mit Statistiken als menschenlesbaren String zurück.
 */
std::string statsAsString(const std::vector<double>& data, const char commentDelimeter[] =
		"\t");

/* I/O */

/*
 * Holt den Input Typgerecht. Von http://www.cplusplus.com/forum/articles/6046/
 * in muss mit dem defaultwert initialisiert sein, bei keiner Angabe wird dieser zurückgegeben.
 */
template<typename T>
T input(T& in) {
	std::string input = "";

	while (true) {
		// füllt den input-string mit der Eingabe
		std::getline(std::cin, input);

		// setzt den Defaultwert bei einem leeren String (getline löscht das letzte \n!)
		if (input == "") {
			return in;
		}

		// Konvertieren von string nach T
		std::stringstream inputStream(input);
		if (inputStream >> in) {
			break;
		}
		std::cout << "Falscher Input! Nochmal." << std::endl;
	}

	return in;
}

/*
 * Konvertiert den gegebenen String in den Typ des gegebenen Arguments. Bei Fehler gibt kommt das Argument zurück.
 */
template<typename T>
T convert(std::string input, T type) {
	// Konvertieren von string nach T
	std::stringstream inputStream(input);
	if (inputStream >> type) {
		return type;
	}

	// Error
	std::cout << "Konnte string '" << input << "' nicht in " << typeid(type).name()
			<< " umwandeln, nehme " << type << "." << std::endl;
	return type;
}

/**
 * Schreibt das Histogramm in eine Datei mit gegebenem Namen, optional mit Statistiken und/oder Header
 *//*
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
}*/

// TODO: write 2D, 3D,

/* Sortiermethoden */

/**
 * Templates um einen Vektor zu (merge-) sortieren. Geändert von http://www.bogotobogo.com/Algorithms/mergesort.php, code A
 */
/**
 *
 */
template<typename T>
std::vector<T> merge(std::vector<T> left, std::vector<T> right) {
	std::vector<T> result;
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
std::vector<T> mergeSort(std::vector<T> m) {
	if (m.size() <= 1)
		return m;

	std::vector<T> left, right, result;
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
// returns eight colors
std::vector<std::string> getColors();

// draw points in geomview
//void drawPointsGeomview(std::vector<std::vector<double> > points, coordinate min = coordinate(-5,-5,-5), coordinate max = coordinate(5,5,5));

/**
 * Plottet ein Histogramm der gegebenen (1D-)Daten mit der Binsize (max-min)/n und der x-Achsen-Beschriftung xlabel.
 */
void plotHist(std::vector<std::vector<double> > data, double min, double max,
		int n, std::vector<std::string> names, const std::string xlabel = "x",
		const std::string file = "/dev/null/");

/**
 * Plottet die Varianz über den Radius. TODO: klären welches vector verwendet wird. ist das von boost besser? müsste immo std sein.
 */
void plotHyperuniformity(std::vector<std::vector<std::vector<double> > > data,
		double xMax, std::vector<std::string> names, const std::string xlabel = "x",
		const std::string ylabel = "y",
		const std::string file = "/dev/null/");

/**
 * Plottet die gegebenen Daten im Format Spalten(Zeilen[3]) mit den gegebenen Grenzen.
 */
void plot3D(std::vector<std::vector<std::vector<double> > > data, std::vector<std::string> names,
		const std::string xlabel = "x", const std::string ylabel = "y",
		const std::string zlabel = "z", const std::string style = "w l");

#endif /* NODELIST_HPP_ */
