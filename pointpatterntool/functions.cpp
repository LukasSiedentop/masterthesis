/*
 * functions.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "functions.hpp"

std::vector<double> stats(std::vector<double> data) {
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
		// check for nan: is not equal to itself
		if (data[i] == data[i]) {
			sum += data[i];
		}
	}

	// expectancy (1st moment) & median
	double expectedValue = sum / data.size();
	double median = data[data.size() / 2];

	// variances (2nd moment)
	double variance = 0;
	double varMedian = 0;

	for (unsigned int i = 0; i < data.size(); i++) {
		// check for nan: is not equal to itself
		if (data[i] == data[i]) {
			variance += pow((data[i] - expectedValue), 2);
			varMedian += pow((data[i] - median), 2);
		}

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

	std::vector<double> statistics;
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

std::string statsAsString(const std::vector<double>& data,
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

	std::stringstream stream;

	// Daten schreiben
	stream << commentDelimeter << "Count: " << data[10] << std::endl;
	stream << commentDelimeter << "Minimum: " << data[11] << ", Maximum: "
			<< data[12] << std::endl;
	stream << commentDelimeter << "Mode: TODO" << std::endl;
	stream << commentDelimeter << "Expectancy +- Standard Deviation: "
			<< data[1] << "+-" << sqrt(data[2]) << ", FWHM: "
			<< 2 * sqrt(2 * log(2) * data[2]) << std::endl;
	stream << commentDelimeter << "Variance: " << data[2] << std::endl;
	stream << commentDelimeter << "Skewness: " << data[3] << " ("
			<< (((data[3] > 0) - (data[3] < 0)) < 0 ?
					"left-skewed" : "right-skewed") << ")" << std::endl;
	stream << commentDelimeter << "Excess Kurtosis: " << data[4] - 3 << " ("
			<< (((data[4] - 3 > 0) - (data[4] - 3 < 0)) < 0 ?
					"platykurtic" : "leptokurtosis ") << ")" << std::endl;
	stream << commentDelimeter << "Median +- Standard Deviation: " << data[13]
			<< "+-" << sqrt(data[14]) << ", FWHM: "
			<< 2 * sqrt(2 * log(2) * data[14]) << std::endl;

	return stream.str();
}

std::vector<std::string> getColors() {
	// colors in HEX #argb
	std::vector<std::string> colors;

	// #55xxxxxx - little bit transparent
	colors.push_back("#55E41A1C"); // red
	colors.push_back("#55377EB8"); // blue
	colors.push_back("#554DAF4A"); // green
	colors.push_back("#55984EA3"); // purple
	colors.push_back("#55FF7F00"); // orange
	colors.push_back("#55FFFF33"); // yellow
	colors.push_back("#55A65628"); // brown
	colors.push_back("#55F781BF"); // pink

	return colors;

}
/*
void drawPointsGeomview(std::vector<std::vector<double> > points, coordinate min, coordinate max) {
	// see http://www.geomview.org/FAQ/answers.shtml
	// initialize geomview for visulaization
	CGAL::Geomview_stream gv(
			CGAL::Bbox_3(min.x(), min.y(), min.z(), max.x(), max.y(), max.z()));
	gv.set_line_width(4);
	// gv.set_trace(true);
	gv.set_bg_color(CGAL::Color(0, 200, 200));
	gv.clear();

	gv << "VECT\n";
	gv << points.size() << " " << points.size() << " 0\n";      // num. of polylines, num. of vertices, num. of colors.

	for (unsigned int i = 0; i < points.size(); i++) {
		gv << "1"; // num. of vertices in each of the 3 polylines,
	}
	gv << "\n";

	for (std::vector<std::vector<double> >::iterator point = points.begin();
			point != points.end(); ++point) {
		gv << (*point)[0] << " " << (*point)[1] << " " << (*point)[2] << "\n"; // Here are the coordinates of each point.
	}

	std::cout << "Continue with Enter." << std::endl;
	std::cin.get();
}*/

void plotHist(std::vector<std::vector<double> > datas, double min, double max, int n,
		std::vector<std::string> names, const std::string xlabel,
		const std::string file) {

	// get filename
	std::stringstream filestream;
	filestream << "tee " << file;
	for (std::vector<std::string>::iterator name = names.begin();
			name != names.end(); ++name) {
		std::string filename = (*name);
		filename.erase(std::remove(filename.begin(), filename.end(), ' '),
				filename.end());
		filestream << "_" << filename;
	}
	filestream << ".gp | gnuplot -persist";

	Gnuplot gp(filestream.str());

	// statistics as comments in file
	for (unsigned i = 0; i < datas.size(); i++) {
		gp << "# " << xlabel << " of pattern '" << names[i] << "'. Statistics: "
				<< std::endl;
		gp << statsAsString(stats(datas[i]), "# ") << std::endl;
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
	gp << "set key left\n";
	gp << "set xlabel '" << xlabel << "'\n";
	gp << "set ylabel 'frequency'\n";

	std::vector<std::string> colors = getColors();
	// build plotstring and send to gnuplot
	std::stringstream plotstring;

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
	for (std::vector<std::vector<double> >::iterator data = datas.begin();
			data != datas.end(); ++data) {
		// The send1d() function sends data to gnuplot's stdin.
		gp.send1d(*data);
	}

	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Continue with Enter." << std::endl;
	std::cin.get();
}

void plotHyperuniformity(std::vector<std::vector<std::vector<double> > > datas, double xMax,
		std::vector<std::string> names, const std::string xlabel,
		const std::string ylabel, const std::string file) {

	double xMin = 0;

	// get filename
	std::stringstream filestream;
	filestream << "tee " << file;
	for (std::vector<std::string>::iterator name = names.begin();
			name != names.end(); ++name) {
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

	std::vector<std::string> colors = getColors();

	std::stringstream plotstring;
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
	for (std::vector<std::vector<std::vector<double> > >::iterator data = datas.begin();
			data != datas.end(); ++data) {

		gp.send1d(*data);
	}

	//gp(fopen("script.gp", "w"));

	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	std::cout << "Continue with Enter." << std::endl;
	std::cin.get();
}

void plot3D(std::vector<std::vector<std::vector<double> > > datas, std::vector<std::string> names,
		const std::string xlabel, const std::string ylabel,
		const std::string zlabel, const std::string style) {
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

	gp << "set key noenhanced\n";

	// Farben
	std::vector<std::string> colors = getColors();

	// Plotstring bauen und an Gnuplot schicken
	std::stringstream plotstring;
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
	for (std::vector<std::vector<std::vector<double> > >::iterator data = datas.begin();
			data != datas.end(); ++data) {

		gp.send1d(*data);
	}

	std::cout << "Continue with Enter." << std::endl;
	std::cin.get();
}

void plot1D(std::vector<std::vector<std::vector<double> > > data, std::vector<std::string> names,
		const std::string xlabel, const std::string ylabel, const std::string style) {
	Gnuplot gp;
		gp << "reset\n";

		gp << "set xlabel '" << xlabel << "'\n";
		gp << "set ylabel '" << ylabel << "'\n";

		gp << "set key noenhanced\n";

		// Farben
		std::vector<std::string> colors = getColors();

		// Plotstring bauen und an Gnuplot schicken
		std::stringstream plotstring;
		plotstring << "plot ";
		for (unsigned i = 0; i < data.size(); i++) {
			plotstring << "'-' " << style << " lt rgb '" << colors[i] << "' t '"
					<< names[i] << "'";
			if (i != data.size() - 1) {
				plotstring << ",";
			}
		}
		gp << plotstring.str() << "\n";

		// Daten senden
		for (std::vector<std::vector<std::vector<double> > >::iterator dat = data.begin();
				dat != data.end(); ++dat) {

			gp.send1d(*dat);
			// TODO: end points ge connected... why?
		}

		std::cout << "Continue with Enter." << std::endl;
		std::cin.get();
}
