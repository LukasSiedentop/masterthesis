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

inline bool coordinate::operator ==(const coordinate &lhs, const coordinate &rhs) {

	bool equals = 0;

	cout << "hier" << endl;
	cout << "jo: " << &this[0] << endl;

	for (unsigned i = 0; i < size(); i++) {
		equals &= abs(rhs[i]-lhs[i]) < tolerance;
	}

	return equals;
}

inline coordinate coordinate::operator +=(const coordinate &lhs, const coordinate &rhs) {

	cout << "hier" << endl;

	for (unsigned i = 0; i < size(); i++) {
		this[i] += rhs[i];
	}

	return *this;
}

