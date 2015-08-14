/*
 * coordinate.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "coordinate.hpp"

using namespace std;

coordinate::coordinate(double point[], unsigned n) {
	for (unsigned i = 0; i < n; i++) {
		push_back(point[i]);
	}
}

coordinate::~coordinate() {
	// TODO: was müssen wir hier tun?
}

int coordinate::numDimensions() {
	return size();
}

bool coordinate::operator ==(const coordinate &rhs) {

	bool equals = 0;

	for (unsigned i = 0; i < size(); i++) {
		equals += abs(rhs[i]-at(i)) < tolerance;
	}

	return equals;
}

bool coordinate::operator !=(const coordinate &rhs) {
	return !(*this == rhs);
}

void coordinate::operator +=(const coordinate &rhs) {
	for (unsigned i = 0; i < size(); i++) {
		at(i) += rhs[i];
	}
}

void coordinate::operator *=(const double &factor) {
	for (unsigned i = 0; i < size(); i++) {
		at(i) *= factor;
	}
}

void coordinate::operator /=(const double &factor) {
	*this *= 1/factor;
}

ostream& operator <<(ostream &os, const coordinate &obj) {
	os << "(";
	for (unsigned i = 0; i < obj.size()-1; i++) {
		os << obj[i] << ", ";
	}
	os << obj.back() << ")";

	return os;
}

double coordinate::lengthSqr() {
	double sumSqr = 0;
	for (unsigned i = 0; i < size(); i++) {
		sumSqr += pow(at(i),2);
	}
	return sumSqr;
}

double coordinate::length() {
	return sqrt(lengthSqr());
}
