/*
 * functions.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "functions.hpp"

vector<double> stats(vector<double> data) {
	// Erwartungswert: E = 1/N \sum_i=0^N x_i
	// Varianz: \sigma^2 = 1/N \sum_i=0^N (x_i - E)^2
	// Schiefe: v = 1/N \sum_i=0^N ((x_i - E)/\sigma)^3
	// Exzess = Wölbung - 3: w =  (1/N \sum_i=0^N ((x_i - E)/\sigma)^4) - 3

	// https://de.wikipedia.org/wiki/Moment_(Stochastik)
	// Moment	Bedeutung
	// 0		=0
	// 1		Erwartungswert
	// 2		Varianz
	// 3		Schiefe
	// 4		Wölbung

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
	double kurtosis = 0;
	for (unsigned int i = 0; i < data.size(); i++) {
		skewness += pow(((data[i] - expectedValue) / sqrt(variance)), 3);
		kurtosis += pow(((data[i] - expectedValue) / sqrt(variance)), 4);
	}
	skewness = skewness / data.size();
	kurtosis = kurtosis / data.size();

	// Vektorposition Bedeutung
	// 0		Modus (TODO) https://de.wikipedia.org/wiki/Modus_(Statistik)
	// 1		Erwartungswert
	// 2		Varianz
	// 3		Schiefe
	// 4		Wölbung
	// 10		Anzahl
	// 11		Minimum
	// 12		Maximum
	// 13		Median
	// 14		variance Median

	vector<double> statistics;
	statistics.resize(15);
	statistics[0] = 0;
	statistics[1] = expectedValue;
	statistics[2] = variance;
	statistics[3] = skewness;
	statistics[4] = kurtosis;
	statistics[10] = data.size();
	statistics[11] = data[0];
	statistics[12] = data[data.size()-1];
	statistics[13] = median;
	statistics[14] = varMedian;

	return statistics;
}

string statsAsString(const vector<double>& data, const char commentDelimeter[]) {
	// https://de.wikipedia.org/wiki/Moment_(Stochastik)
	// Moment	Bedeutung
	// 0		=0
	// 1		Erwartungswert
	// 2		Varianz
	// 3		Schiefe
	// 4		Wölbung

	// Vektorposition Bedeutung
	// 0		Modus (TODO) https://de.wikipedia.org/wiki/Modus_(Statistik)
	// 1		Erwartungswert
	// 2		Varianz
	// 3		Schiefe
	// 4		Wölbung
	// 10		Anzahl
	// 11		Minimum
	// 12		Maximum
	// 13		Median
	// 14		variance Median

	stringstream stream;

	// Daten schreiben
	stream << commentDelimeter << "Anzahl: " << data[10] << endl;
	stream << commentDelimeter << "Minimum: " << data[11] << ", Maximum: "
			<< data[12] << endl;
	stream << commentDelimeter << "Modus: TODO" << endl;
	stream << commentDelimeter << "Erwartungswert +- Standardabweichung: "
			<< data[1] << "+-" << sqrt(data[2]) << ", FWHM: "
			<< 2 * sqrt(2 * log(2) * data[2]) << endl;
	stream << commentDelimeter << "Varianz: " << data[2] << endl;
	stream << commentDelimeter << "Schiefe: " << data[3] << " ("
			<< (((data[3] > 0) - (data[3] < 0)) < 0 ?
					"linksschief" : "rechtschief") << ")" << endl;
	stream << commentDelimeter << "Exzess: " << data[4] - 3 << " ("
			<< (((data[4] - 3 > 0) - (data[4] - 3 < 0)) < 0 ?
					"flachgipflig" : "steilgipflig") << ")" << endl;
	stream << commentDelimeter << "Median +- Standardabweichung: " << data[13]
			<< "+-" << sqrt(data[14]) << ", FWHM: "
			<< 2 * sqrt(2 * log(2) * data[14]) << endl;

	return stream.str();
}

void plot3D(vector<vector<coordinate> > datas, const char xlabel[],
		const char ylabel[], const char zlabel[], const char style[]) {
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
		plotstring << "'-' " << style << " lt rgb '" << colors[i] << "' t 'Muster " << i << "'";
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
