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
#include <limits>
#include <cstdlib>
#include <iterator>
#include <stdio.h>
/*
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K>      Triangulation;
*/
static const double tolerance = 0.0000001;

/**
 * Represents a point in n-Dimensional space with various operators
 */
class coordinate {
private:
	std::vector<double> position;
public:
	// Reference: https://en.wikipedia.org/wiki/Rule_of_three_(C%2B%2B_programming)
	coordinate();
	// Generates a random point in n dimensions with all coordinates the same value
	coordinate(double pos, unsigned int n);
	coordinate(double x, double y);
	coordinate(double x, double y, double z);
	// Generates a random point with coordinate in (0,1)^n
	coordinate(unsigned int n);

	std::vector<double>* getVector();

	unsigned int dimensions() const;

	// subscript operator
	const double& operator[](const int i) const;
	double& operator[](const int i);

	// Comparision: all components need to be equal
	bool operator ==(const coordinate& rhs);
	bool operator !=(const coordinate& rhs);

	// Gößer/Kleiner über Länge?
	//bool operator <(const coordinate &rhs);
	//bool operator >(const coordinate &rhs);
	//bool operator <=(const coordinate &rhs);
	//bool operator >=(const coordinate &rhs);

	// Shift this point with rhs.
	coordinate & operator +=(const coordinate& rhs);
	// Shift this point with summand in every component.
	coordinate & operator +=(const double& summand);
	// Shift this point with -rhs.
	coordinate & operator -=(const coordinate& rhs);
	// Shift this point with -subtrahend in every component.
	coordinate & operator -=(const double& subtrahend);
	// Multiply this point componentwise with rhs TODO: necessary? implement dot-product instead
	coordinate & operator *=(const coordinate& rhs);
	// Scale this point with factor
	coordinate & operator *=(const double& factor);
	// Scale this point with 1/factor
	coordinate & operator /=(const double& factor);

	// Gives an adjustable string-representation of the point
	std::string toString(const std::string begin = "(", const std::string delimiter = ",",
			const std::string end = ")", const int precision = 4) const;
	// Representation of the point in a stringstream
	friend std::ostream& operator <<(std::ostream& os, const coordinate& obj);

	// Adds lhs's to rhs's components.
	friend coordinate operator+(const coordinate& lhs, const coordinate& rhs);
	// Adds a summand to every component.
	friend coordinate operator+(const coordinate& lhs, const double& summand);
	// Subtracts lhs's components from rhs's.
	friend coordinate operator-(const coordinate& lhs, const coordinate& rhs);
	// Subtracts a subtrahend from every component.
	friend coordinate operator-(const coordinate& lhs,
			const double& subtrahend);
	// Multiplies the components of both points.
	friend coordinate operator*(const coordinate& lhs, const coordinate& rhs);
	// Multiplies each component with a factor.
	friend coordinate operator*(const coordinate& lhs, const double& factor);
	// Divides each component by a divisor.
	friend coordinate operator/(const coordinate& lhs, const double& divisor);

	// checks whether this coordinate resides within the given axis aligned bounding box
	bool insideAABB(const coordinate& min, const coordinate& max) const;
	// checks whether this point is within the given cube
	bool insideAABB(const coordinate& mid, const double& r) const;

	// get positions
	double x() const;
	double y() const;
	double z() const;

	// Returns the minimal component of the coordinate.
	double min() const;
	// Returns the maximal component of the coordinate.
	double max() const;

	double euklidian(const coordinate& point) const;
	double lengthSqr() const;
	double length() const;

	// returns the angle between this and the given vector in radians
	double angle(const coordinate& b) const;

	// returns the shifters from this coordinate
	std::vector<coordinate> getShifters();

	// returns the CGAL point
	//Triangulation::Point getPoint() const;

	// dot-product etween two vectors this and b
	//static double scp(const coordinate& a, const coordinate& b);
	double scp(const coordinate& b) const;

	// cross-product
	static coordinate cpr(const coordinate& a, const coordinate& b);

	// TODO: unnecessary when stored shifter vector is used
	// Gives the shortest vector from a to b, considering the given shifting vectors for periodic boundary conditions.
	//static coordinate getVec(const coordinate& a, const coordinate& b, std::vector<coordinate> shifters);
};

#endif /* COORDINATE_HPP_ */
