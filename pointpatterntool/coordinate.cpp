/*
 * coordinate.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "coordinate.hpp"

using namespace std;

coordinate::coordinate(double x, double y) {
	push_back(x);
	push_back(y);
}

coordinate::coordinate(double x, double y, double z) {
	push_back(x);
	push_back(y);
	push_back(z);
}

coordinate::coordinate(double point[], unsigned n) {
	for (unsigned i = 0; i < n; i++) {
		push_back(point[i]);
	}
}

coordinate::~coordinate() {
	// TODO: was müssen wir hier tun? Nötig?
	vector<double>().swap(*this);
}

int coordinate::numDimensions() {
	return size();
}

coordinate coordinate::operator+(const coordinate & lhs,const coordinate & rhs) {
	return coordinate(lhs) += rhs;
}

coordinate coordinate::operator-(const coordinate & lhs, const coordinate & rhs) {
	return coordinate(lhs) -= rhs;
}

coordinate coordinate::operator*(const coordinate & lhs, const coordinate & rhs) {
	return coordinate(lhs) *= rhs;
}

coordinate coordinate::operator/(const coordinate & lhs, const coordinate & rhs) {
	return coordinate(lhs) /= rhs;
}

bool coordinate::operator ==(const coordinate &rhs) {

	bool equals = 0;

	//TODO vergleich mit
	//for (vector<double>::iterator it = begin(); it != end(); ++it) {}
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

void coordinate::operator -=(const coordinate &rhs) {
	for (unsigned i = 0; i < size(); i++) {
		at(i) -= rhs[i];
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

std::string coordinate::toString(const char begin[]="(", const char delimiter[]=", ", const char end[]=")") {
	stringstream stream;
	stream << begin;

	for (unsigned i = 0; i < size()-1; i++) {
		stream << at(i) << delimiter;
	}

	stream << back() << end;

	return stream.str();
}

double coordinate::euklidian(coordinate & point) {
	// TODO: geht vielleicht nicht...
	return coordinate(this-point).length();
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

static double coordinate::scp(const coordinate &a, const coordinate &b) {
	double scp = 0;
	for (unsigned i = 0; i < a.size(); i++) {
		scp += a[i]*b[i];
	}
	return scp;
}


static coordinate coordinate::getVec(const coordinate & a, const coordinate & b, vector<coordinate> shifters) {
	coordinate vec = a-b;
	double lengthSqr = vec.lengthSqr();

	// verschobener Vektor
	coordinate tmpVec;

	for (unsigned i = 0; i < shifters.size(); i++) {
		tmpVec = vec+shifters[i];

		// vergleich der Länge mit dem verschobenen Vektor
		if(tmpVec.lengthSqr() < lengthSqr) {
			vec += shifters[i];
		}
	}

	return vec;
}
