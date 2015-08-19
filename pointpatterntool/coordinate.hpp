/*
 * coordinate.hpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#ifndef COORDINATE_HPP_
#define COORDINATE_HPP_

#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>

#include <cstdlib>

/**
 * Point in n-Dimensional space with various operators
 */
class coordinate: public std::vector<double> {
private:
	static const double tolerance = 0.0000001;
public:
	/* Standardkonstruktor */
	coordinate();
	/* Züfällige Koordinate in der Box (0,1)^n*/
	coordinate(unsigned int n);
	/* Copy Constructor*/
	coordinate(const coordinate & coord);
	/* 2D Punkt aus Koordinaten */
	coordinate(double x, double y);
	/* 3D Punkt aus Koordinaten */
	coordinate(double x, double y, double z);
	/* nD Punkt aus Koordinaten */
	coordinate(double point[], unsigned n);
	/* Destruktor */
	~coordinate();

	/* Anzahl der Dimensionen */
	int numDimensions();

	/* Vergleichsoperatoren: alle Komponenten müssen gleich sein */
	bool operator ==(const coordinate &rhs);
	/* negierter Vergleichsoperatoren: mindestens eine Komponente ist ungleich */
	bool operator !=(const coordinate &rhs);
	// Gößer/Kleiner über Länge?
	//bool operator <(const coordinate &rhs);
	//bool operator >(const coordinate &rhs);
	//bool operator <=(const coordinate &rhs);
	//bool operator >=(const coordinate &rhs);

	/* verschiebt diese Koordinate um rhs */
	coordinate & operator +=(const coordinate &rhs);
	/* verschiebt diese Koordinate um -rhs */
	coordinate & operator -=(const coordinate &rhs);
	/* multipliziert diese Koordinate Komponentenweise mit dem gegebenen */
	coordinate & operator *=(const coordinate &rhs);
	/* skaliert diesen Vektor um einen factor */
	coordinate & operator *=(const double &factor);
	/* skaliert diesen Vektor um einen 1/factor */
	coordinate & operator /=(const double &factor);

	/* Gibt eine String-Repräsentation der Koordinate mit anpassbaren Trennzeichen */
	std::string toString(const char begin[] = "(", const char delimiter[] = ",",
			const char end[] = ")");
	/* Repräsentation der Koordinate im Stringstream */
	friend std::ostream& operator <<(std::ostream &os, const coordinate &obj);

	/* Addiert die Vektoren Komponentenweise. Teuer! */
	friend coordinate operator+(const coordinate & lhs, const coordinate & rhs);
	/* Subtrahiert die Vektoren Komponentenweise. Teuer! */
	friend coordinate operator-(const coordinate & lhs, const coordinate & rhs);
	/* multipliziert jede Komponente mit factor */
	friend coordinate operator*(const coordinate & lhs, const coordinate & rhs);
	/* multipliziert jede Komponente mit factor */
	friend coordinate operator*(const coordinate & lhs, const double &factor);
	/* dividiert jede Komponente mit factor */
	friend coordinate operator/(const coordinate & lhs, const double &factor);

	/* Gibt den euklidschen Abstand zwischen dieser Koordinate und point zurück */
	double euklidian(coordinate point);
	/* Gibt die quadrierte euklidsche Länge des Vektors vom Ursprung zu dieser Koordinate zurück. Etwa doppelt so schnell wie length() */
	double lengthSqr();
	/* Gibt die euklidsche Länge des Vektors vom Ursprung zu dieser Koordinate zurück */
	double length();

	/* Gibt das Skalarprodukt der gegebenen Koordinaten zurück */
	static double scp(const coordinate &a, const coordinate &b);

	/* Gibt den minimal langen Vektor von A nach B zurück, unter Berücksichtigung periodischer Randbedingungen */
	static coordinate getVec(const coordinate & a, const coordinate & b,
			vector<coordinate> shifters);
};

#endif /* COORDINATE_HPP_ */
