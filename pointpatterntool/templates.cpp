/*
 * templates.cpp
 *
 *  Created on: 03.08.2015
 *      Author: lukas
 */

#include <vector>

#include "gnuplot-iostream.h"

using namespace std;

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


template<typename T>
void plotHist(vector<T> data, double min, double max, int n, const char xlabel[]){
	Gnuplot gp;
	// Don't forget to put "\n" at the end of each line!

	gp << "reset\n";

	gp << "n=" << n << "\n";
	gp << "max=" << max << "\n";
	gp << "min=" << min << "\n";
	gp << "width=" << (max-min)/n << "\n";

	gp << "hist(x,width)=width*floor(x/width) #+width/2.0\n";

	gp << "set xrange [min:max]\n";
	gp << "set yrange [0:]\n";
	gp << "set offset graph 0.05,0.05,0.05,0.0\n";
	// 10 x-tics
	gp << "set xtics min," << (max-min)/10 << ",max\n";
	gp << "set boxwidth width*0.9\n";
	gp << "set style fill solid 0.5\n";
	gp << "set tics out nomirror\n";
	gp << "set xlabel '" << xlabel << "'\n";
	gp << "set ylabel 'Häufigkeit'\n";

	// '-' means read from stdin.  The send1d() function sends data to gnuplot's stdin.
	gp << "plot '-' u ($1):(width) w boxes smooth freq lc rgb'blue' notitle\n";
	gp.send1d(data);

	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Weiter mit Enter." << std::endl;
	std::cin.get();
}
