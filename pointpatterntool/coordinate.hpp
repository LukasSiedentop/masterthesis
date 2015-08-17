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

/**
 * Point in n-Dimensional space with various operators
 */
class coordinate: public std::vector<double> {
private:
	static const double tolerance = 2e-6;
public:
	coordinate();
	coordinate(const coordinate & coord);
	coordinate(double x, double y);
	coordinate(double x, double y, double z);
	coordinate(double point[], unsigned n);
	~coordinate();

	/* Anzahl der Dimensionen */
	int numDimensions();

	/* boolsche Operatoren */
	bool operator ==(const coordinate &rhs);
	bool operator !=(const coordinate &rhs);
	// Gößer/Kleiner über Länge?
	//bool operator <(const coordinate &rhs);
	//bool operator >(const coordinate &rhs);
	//bool operator <=(const coordinate &rhs);
	//bool operator >=(const coordinate &rhs);

	/* bewegen + skalieren */
	coordinate & operator +=(const coordinate &rhs);
	coordinate & operator -=(const coordinate &rhs);
	coordinate & operator *=(const double &factor);
	coordinate & operator /=(const double &factor);

	/* ausgeben */
	std::string toString(const char begin[] = "(", const char delimiter[] = ",",
			const char end[] = ")");
	friend std::ostream& operator <<(std::ostream &os, const coordinate &obj);

	/* Addition/Subtraktion/Multiplikation/Division */
	friend coordinate operator+(const coordinate & lhs, const coordinate & rhs);
	friend coordinate operator-(const coordinate & lhs, const coordinate & rhs);
	friend coordinate operator*(const coordinate & lhs, const double &rhs);
	friend coordinate operator/(const coordinate & lhs, const double &rhs);

	/* abstände */
	double euklidian(coordinate & point);
	double lengthSqr();
	double length();

	/* skalarprodukt */
	static double scp(const coordinate &a, const coordinate &b);

	/* gibt den minimal Langen Vektor von A nach B zurück, unter Berücksichtigung periodischer Randbedingungen */
	static coordinate getVec(const coordinate & a, const coordinate & b,
			vector<coordinate>  shifters);
};

#endif /* COORDINATE_HPP_ */
