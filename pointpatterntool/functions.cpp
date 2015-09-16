/*
 * functions.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "functions.hpp"

vector<double> stats(vector<double> data) {
	// expectancy: E = 1/N \sum_i=0^N x_i
	// variance: \sigma^2 = 1/N \sum_i=0^N (x_i - E)^2
	// skewness: v = 1/N \sum_i=0^N ((x_i - E)/\sigma)^3
	// kurtosis excess = kurtosis - 3: w =  (1/N \sum_i=0^N ((x_i - E)/\sigma)^4) - 3

	// https://de.wikipedia.org/wiki/Moment_(Stochastik)
	// moment	meaning
	// 0		=0
	// 1		expectancy
	// 2		variance
	// 3		skewness
	// 4		kurtosis

	// sum
	double sum = 0;
	for (unsigned int i = 0; i < data.size(); i++) {
		sum += data[i];
	}

	// expectancy (1st moment) & median
	double expectedValue = sum / data.size();
	double median = data[data.size() / 2];

	// variances (2nd moment)
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

	/*
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

	vector<double> statistics;
	statistics.resize(15);
	statistics[0] = 0;
	statistics[1] = expectedValue;
	statistics[2] = variance;
	statistics[3] = skewness;
	statistics[4] = kurtosis;
	statistics[10] = data.size();
	statistics[11] = data[0];
	statistics[12] = data[data.size() - 1];
	statistics[13] = median;
	statistics[14] = varMedian;

	return statistics;
}

string statsAsString(const vector<double>& data,
		const char commentDelimeter[]) {
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
	stream << commentDelimeter << "Count: " << data[10] << endl;
	stream << commentDelimeter << "Minimum: " << data[11] << ", Maximum: "
			<< data[12] << endl;
	stream << commentDelimeter << "Mode: TODO" << endl;
	stream << commentDelimeter << "Expectancy +- Standard Deviation: "
			<< data[1] << "+-" << sqrt(data[2]) << ", FWHM: "
			<< 2 * sqrt(2 * log(2) * data[2]) << endl;
	stream << commentDelimeter << "Variance: " << data[2] << endl;
	stream << commentDelimeter << "Skewness: " << data[3] << " ("
			<< (((data[3] > 0) - (data[3] < 0)) < 0 ?
					"left-skewed" : "right-skewed") << ")" << endl;
	stream << commentDelimeter << "Excess Kurtosis: " << data[4] - 3 << " ("
			<< (((data[4] - 3 > 0) - (data[4] - 3 < 0)) < 0 ?
					"platykurtic" : "leptokurtosis ") << ")" << endl;
	stream << commentDelimeter << "Median +- Standard Deviation: " << data[13]
			<< "+-" << sqrt(data[14]) << ", FWHM: "
			<< 2 * sqrt(2 * log(2) * data[14]) << endl;

	return stream.str();
}

vector<string> getColors() {
	// colors in HEX #argb
	vector<string> colors;

	colors.push_back("#55E41A1C");
	colors.push_back("#55377EB8");
	colors.push_back("#554DAF4A");
	colors.push_back("#55984EA3");
	colors.push_back("#55FF7F00");
	colors.push_back("#55FFFF33");
	colors.push_back("#55A65628");
	colors.push_back("#55F781BF");

	return colors;

}

void plotHist(vector<vector<double> > datas, double min, double max, int n,
		vector<string> names, const char xlabel[], const char file[]) {

	// get filename
	stringstream filestream;
	filestream << "tee " << file;
	for (vector<string>::iterator name = names.begin(); name != names.end();
			++name) {
		string filename = (*name);
		filename.erase(std::remove(filename.begin(), filename.end(), ' '),
				filename.end());
		filestream << "_" << filename;
	}
	filestream << ".gp | gnuplot -persist";

	Gnuplot gp(filestream.str());

	// statistics as comments in file
	for (unsigned i = 0; i < datas.size(); i++) {
		gp << "# " << xlabel << " of pattern '" << names[i] << "'. Statistics: "
				<< endl;
		gp << statsAsString(stats(datas[i]), "# ") << endl;
	}

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
	gp << "set ylabel 'frequency'\n";

	vector<string> colors = getColors();
	// build plotstring and send to gnuplot
	stringstream plotstring;

	plotstring << "plot ";
	for (unsigned i = 0; i < datas.size(); i++) {

		// '-' means read from stdin
		plotstring << "'-' u (hist($1, width)):(1.0/" << datas[i].size()
				<< ") w boxes smooth freq lc rgb '" << colors[i] << "' t '"
				<< names[i] << "'";
		if (i != datas.size() - 1) {
			plotstring << ",";
		}
	}
	gp << plotstring.str() << "\n";

	// Daten senden
	for (vector<vector<double> >::iterator data = datas.begin();
			data != datas.end(); ++data) {
		// The send1d() function sends data to gnuplot's stdin.
		gp.send1d(*data);
	}

	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	cout << "Weiter mit Enter." << endl;
	cin.get();
}

void plotHyperuniformity(vector<vector<vector<double> > > datas, double xMax,
		vector<string> names, const char xlabel[], const char ylabel[],
		string file) {

	double xMin = 0;

	// get filename
	stringstream filestream;
	filestream << "tee " << file;
	for (vector<string>::iterator name = names.begin(); name != names.end();
			++name) {
		std::cout << (*name) << "\n";
		(*name).erase(std::remove((*name).begin(), (*name).end(), ' '),
				(*name).end());
		std::cout << (*name) << "\n";
		filestream << "_" << *name;
	}
	filestream << ".gp | gnuplot -persist";

	Gnuplot gp(filestream.str());
	gp << "reset\n";

	gp << "set fit logfile '/dev/null'\n";

	gp << "set key top left\n";
	gp << "set xrange [" << xMin << ":" << xMax << "]\n";
	gp << "set yrange [0:]\n";
	//gp << "set offset graph 0.05,0.05,0.05,0.0\n";
	// 10 x-tics
	gp << "set xtics " << xMin << ",1," << xMax << "\n";
	gp << "set tics out nomirror\n";
	gp << "set xlabel '" << xlabel << "'\n";
	gp << "set ylabel '" << ylabel << "'\n";

	vector<string> colors = getColors();

	stringstream plotstring;
	for (unsigned i = 0; i < datas.size(); i++) {

		gp << "A" << i << " = 1\n";
		gp << "B" << i << " = 2\n";
		gp << "f" << i << "(x) = A" << i << "*x**B" << i << "\n";
		gp << "fit [0:" << xMax << "] f" << i << "(x) '-' u 1:2 via A" << i
				<< "\n"; // , B" << i << "
		gp.send1d(datas[i]);

		if (i == 0) {
			plotstring << "plot ";
		}

		// '-' means read from stdin.  The send1d() function sends data to gnuplot's stdin.
		plotstring << "'-' u 1:2 ls 7 lc rgb'" << colors[i] << "' t 'variance "
				<< names[i] << "', f" << i << "(x) lc rgb'" << colors[i]
				<< "' t sprintf('fit: f_" << i << "(R) = %.3fR^{%.1f}',A" << i // {\264}
				<< ", B" << i << ")";

		if (i != datas.size() - 1) {
			plotstring << ",";
		}
	}

	gp << plotstring.str() << "\n";

	// Daten senden
	for (vector<vector<vector<double> > >::iterator data = datas.begin();
			data != datas.end(); ++data) {

		gp.send1d(*data);
	}

	//gp(fopen("script.gp", "w"));

	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Weiter mit Enter." << std::endl;
	std::cin.get();
}

void plot3D(vector<vector<vector<double> > > datas, vector<string> names,
		const char xlabel[], const char ylabel[], const char zlabel[],
		const char style[]) {
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
	vector<string> colors = getColors();

	// Plotstring bauen und an Gnuplot schicken
	stringstream plotstring;
	plotstring << "splot ";
	for (unsigned i = 0; i < datas.size(); i++) {
		plotstring << "'-' " << style << " lt rgb '" << colors[i] << "' t '"
				<< names[i] << "'";
		if (i != datas.size() - 1) {
			plotstring << ",";
		}
	}
	gp << plotstring.str() << "\n";

	// Daten senden
	for (vector<vector<vector<double> > >::iterator data = datas.begin();
			data != datas.end(); ++data) {

		gp.send1d(*data);
	}

	cout << "Weiter mit Enter." << endl;
	cin.get();
}
