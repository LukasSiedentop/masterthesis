/*
 * coordinate.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "coordinate.hpp"

using namespace std;
coordinate::coordinate() {
}

coordinate::coordinate(unsigned int n) {
	for (unsigned int i = 0; i < n; i++) {
		push_back((double) rand() / RAND_MAX);
	}
}

coordinate::coordinate(const coordinate& coord) {
	*this = coord;
}

coordinate::coordinate(double x, double y) {
	push_back(x);
	push_back(y);
}

coordinate::coordinate(double x, double y, double z) {
	/*this->resize(3);
	this[0]=x;
	this[1]=y;
	this[2]=z;*/

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

bool coordinate::operator ==(const coordinate &rhs) {
	for (unsigned i = 0; i < size(); i++) {
		if (abs(rhs[i] - (*this)[i]) > tolerance) {
			return false;
		}
	}
	return true;
}

bool coordinate::operator !=(const coordinate &rhs) {
	return !(*this == rhs);
}

coordinate & coordinate::operator +=(const coordinate &rhs) {
	for (unsigned i = 0; i < size(); i++) {
		(*this)[i] += rhs[i];
	}

	return *this;
}

coordinate & coordinate::operator +=(const double &summand) {
	for (unsigned i = 0; i < size(); i++) {
		(*this)[i] -= summand;
	}
	return *this;
}

coordinate & coordinate::operator -=(const coordinate &rhs) {
	for (unsigned i = 0; i < size(); i++) {
		(*this)[i] -= rhs[i];
	}
	return *this;
}

coordinate & coordinate::operator -=(const double &subtrahend) {
	for (unsigned i = 0; i < size(); i++) {
		(*this)[i] -= subtrahend;
	}
	return *this;
}

coordinate & coordinate::operator *=(const coordinate &rhs) {
	for (unsigned i = 0; i < size(); i++) {
		(*this)[i] *= rhs[i];
	}
	return *this;
}

coordinate & coordinate::operator *=(const double &factor) {
	for (unsigned i = 0; i < size(); i++) {
		(*this)[i] *= factor;
	}
	return *this;
}

coordinate & coordinate::operator /=(const double &factor) {
	*this *= 1 / factor;
	return *this;
}

ostream& operator <<(ostream &os, const coordinate &obj) {
	os << "(";
	for (unsigned i = 0; i < obj.size() - 1; i++) {
		os << obj[i] << ", ";
	}
	os << obj.back() << ")";

	return os;
}

string coordinate::toString(const char begin[], const char delimiter[],
		const char end[]) {
	stringstream stream;
	stream << begin;

	for (unsigned i = 0; i < size() - 1; i++) {
		stream << (*this)[i] << delimiter;
	}

	stream << back() << end;

	return stream.str();
}

double coordinate::min() {
	double minimu = numeric_limits<double>::infinity();
	for (unsigned i = 0; i < size(); i++) {
		minimu = std::min(minimu, (*this)[i]);
	}
	return minimu;
}

double coordinate::euklidian(coordinate point) {
	return coordinate(*this - point).length();
}

double coordinate::lengthSqr() const {
	double sumSqr = 0;
	for (unsigned i = 0; i < size(); i++) {
		sumSqr += (*this)[i] * (*this)[i];
	}
	return sumSqr;
}

double coordinate::length() {
	return sqrt(lengthSqr());
}

double coordinate::scp(const coordinate &a, const coordinate &b) {
	double scp = 0;
	for (unsigned i = 0; i < a.size(); i++) {
		scp += a[i] * b[i];
	}
	return scp;
}

coordinate coordinate::getVec(const coordinate & a, const coordinate & b,
		vector<coordinate> shifters) {
	coordinate vec = a - b;
	double maxFeatureSize = shifters.back().lengthSqr() / 4;

	// verschobener Vektor
	coordinate tmpVec;

	for (unsigned i = 0; i < shifters.size(); i++) {
		// Vektor verschieben
		tmpVec = vec + shifters[i];

		// vergleich der Länge mit dem verschobenen Vektor
		if (tmpVec.lengthSqr() < maxFeatureSize) {
			return tmpVec;
		}
	}

	// Error
	return vec;
}

coordinate operator+(const coordinate& lhs, const coordinate& rhs) {
	return coordinate(lhs) += rhs;
}

coordinate operator+(const coordinate& lhs, const double& summand) {
	return coordinate(lhs) += summand;
}

coordinate operator-(const coordinate& lhs, const coordinate& rhs) {
	return coordinate(lhs) -= rhs;
}

coordinate operator-(const coordinate& lhs, const double& subtrahend) {
	return coordinate(lhs) -= subtrahend;
}

coordinate operator*(const coordinate& lhs, const coordinate& rhs) {
	return coordinate(lhs) *= rhs;
}

coordinate operator*(const coordinate& lhs, const double& factor) {
	return coordinate(lhs) *= factor;
}

coordinate operator/(const coordinate& lhs, const double& factor) {
	return coordinate(lhs) /= factor;
}
