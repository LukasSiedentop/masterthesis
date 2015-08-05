/*
 * templates.cpp
 *
 *  Created on: 03.08.2015
 *      Author: lukas
 */

#include <vector>

#include <string>
#include <sstream>

#include "gnuplot-iostream.h"

using namespace std;

/*
 * Holt den Input Typgerecht. Von http://www.cplusplus.com/forum/articles/6046/
 * In muss mit dem defaultwert initialisiert sein, bei keiner Angabe wird dieser zurückgegeben.
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

/**
 * Berechnet und schreibt die Statistiken eines Vektors auf die Konsole. Dieser muss aus Daten bestehen, die Darstell-, Addier-, und Teilbar sind.
 * Zurzeit: Anzahl, Erwartungswert mit Standardabweichung und FWHM, Median mit Standardabweichung und FWHM
 */
template<typename T>
void stats(vector<T> data) {
	// Summe
	double sum = 0;
	for (unsigned int i = 0; i < data.size(); i++) {
		sum += data[i];
	}

	// Erwartungswert (1. Moment) & Median
	double expectedValue = sum / data.size();
	double median = data[data.size() / 2];

	// Varianzen (2. Moment)
	double variance = 0;
	double varMedian = 0;

	for (unsigned int i = 0; i < data.size(); i++) {
		variance += pow((data[i] - expectedValue), 2);
		varMedian += pow((data[i] - median), 2);
	}
	variance = variance / data.size();
	varMedian = varMedian / data.size();

	double skewness = 0;
	for (unsigned int i = 0; i < data.size(); i++) {
		skewness += pow(((data[i] - expectedValue) / sqrt(variance)), 3);
	}
	skewness = skewness / data.size();

	// https://de.wikipedia.org/wiki/Moment_(Stochastik)
	// Moment	Bedeutung
	// 0		=0
	// 1		Erwartungswert
	// 2		Varianz
	// 3		Schiefe (<0: linksschief, >0: rechtsschief)
	// 4		Wölbung

	// TODO: Modus

	// Daten schreiben
	cout << char(9) << "Anzahl: " << data.size() << endl;
	cout << char(9) << "Minimum: " << data[0] << ", Maximum: "
			<< data[data.size() - 1] << endl;
	cout << char(9) << "Erwartungswert +- Standardabweichung: " << expectedValue
			<< "+-" << sqrt(variance) << ", FWHM: "
			<< 2 * sqrt(2 * log(2) * variance) << endl;
	cout << char(9) << "Varianz: " << variance << endl;
	cout << char(9) << "Schiefe: " << skewness << endl;
	cout << char(9) << "Median +- Standardabweichung: " << median << "+-"
			<< sqrt(varMedian) << ", FWHM: " << 2 * sqrt(2 * log(2) * varMedian)
			<< endl;
}

/**
 * Templates um einen Vektor zu (merge-) sortieren. Geändert von http://www.bogotobogo.com/Algorithms/mergesort.php, code A
 * TODO: umbaubar auf die NodeList? Diese als Vector?
 */
template<typename T>
vector<T> merge(vector<T> left, vector<T> right) {
	vector<T> result;
	while ((int) left.size() > 0 || (int) right.size() > 0) {
		if ((int) left.size() > 0 && (int) right.size() > 0) {
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

/**
 * Plottet ein Histogramm der gegebenen (1D-)Daten mit der Binsize (max-min)/n und der x-Achsen-Beschriftung xlabel.
 */
template<typename T>
void plotHist(vector<T> data, double min, double max, int n,
		const char xlabel[]) {
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
void plot2D(vector<vector<T> > data) {

	Gnuplot gp;
	gp << "reset\n";

	gp << "A = 1\n";
	gp << "f(x) = A*x**2\n";
	gp << "fit [0:10] f(x) '-' u 1:2 via A\n";
	gp.send1d(data);

	gp << "set xrange [0:10]\n";
	gp << "set yrange [0:]\n";
	//gp << "set offset graph 0.05,0.05,0.05,0.0\n";
	// 10 x-tics
	gp << "set xtics 0,1,10\n";
	gp << "set tics out nomirror\n";
	gp << "set xlabel 'Radius R'\n";
	gp << "set ylabel 'Varianz  {/Symbol s}^2(R)'\n";

	// '-' means read from stdin.  The send1d() function sends data to gnuplot's stdin.
	gp
			<< "plot '-' u 1:2 ls 7 lc rgb'blue' t 'Numerische Daten', f(x) lc rgb'blue' t sprintf('Fit: f(R) = %dR^2',A)\n";
	gp.send1d(data);

	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Weiter mit Enter." << std::endl;
	std::cin.get();
}
