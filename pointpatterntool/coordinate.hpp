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

/**
 * Point in n-Dimensional space with various operators
 */
class coordinate: std::vector<double> {
private:
	static const double tolerance = 2e-6;
public:
	coordinate();
	coordinate(double point[], unsigned n);
	~coordinate();

	/* nützliche Operatoren */
	bool operator ==(const coordinate &rhs);
	bool operator !=(const coordinate &rhs);
	void operator +=(const coordinate &rhs);
	void operator *=(const double &factor);
	void operator /=(const double &factor);
	friend std::ostream& operator <<(std::ostream &os, const coordinate &obj);

	double lengthSqr();
	double length();

	int numDimensions();
};

#endif /* COORDINATE_HPP_ */
