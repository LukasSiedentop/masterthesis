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
	position.resize(n);
	for (unsigned int i = 0; i < n; i++) {
		position[i] = (double) rand() / RAND_MAX;
	}
}

coordinate::coordinate(double x, double y) {
	position.resize(2);
	position[0] = x;
	position[1] = y;
}

coordinate::coordinate(double x, double y, double z) {
	position.resize(3);
	position[0] = x;
	position[1] = y;
	position[2] = z;
}

vector<double>* coordinate::getVector() {
	return &position;
}

unsigned int coordinate::dimensions() const {
	return position.size();
}

double& coordinate::operator [](const int i){
	return position[i];
}

const double& coordinate::operator [](const int i) const{
	return position[i];
}

bool coordinate::operator ==(const coordinate &rhs) {
	for (unsigned int i = 0; i < dimensions(); i++) {
		if (std::abs(rhs[i] - position[i]) > tolerance) {
			return false;
		}
	}
	return true;
}

bool coordinate::operator !=(const coordinate &rhs) {
	return !(*this == rhs);
}

coordinate & coordinate::operator +=(const coordinate &rhs) {
	for (unsigned i = 0; i < dimensions(); i++) {
		(*this)[i] += rhs[i];
	}

	return *this;
}

coordinate & coordinate::operator +=(const double &summand) {
	for (unsigned i = 0; i < dimensions(); i++) {
		(*this)[i] -= summand;
	}
	return *this;
}

coordinate & coordinate::operator -=(const coordinate &rhs) {
	for (unsigned i = 0; i < dimensions(); i++) {
		(*this)[i] -= rhs[i];
	}
	return *this;
}

coordinate & coordinate::operator -=(const double &subtrahend) {
	for (unsigned i = 0; i < dimensions(); i++) {
		(*this)[i] -= subtrahend;
	}
	return *this;
}

coordinate & coordinate::operator *=(const coordinate &rhs) {
	for (unsigned i = 0; i < dimensions(); i++) {
		(*this)[i] *= rhs[i];
	}
	return *this;
}

coordinate & coordinate::operator *=(const double &factor) {
	for (unsigned i = 0; i < dimensions(); i++) {
		(*this)[i] *= factor;
	}
	return *this;
}

coordinate & coordinate::operator /=(const double &factor) {
	*this *= 1 / factor;
	return *this;
}

ostream& operator <<(ostream &os, const coordinate &obj) {
	return os << obj.toString();
}

string coordinate::toString(const char begin[], const char delimiter[],
		const char end[]) const {
	stringstream stream;
	stream << begin;

	for (unsigned i = 0; i < dimensions() - 1; i++) {
		stream << (*this)[i] << delimiter;
	}

	stream << position[dimensions()-1] << end;

	return stream.str();
}

double coordinate::min() const{
	double minimum = numeric_limits<double>::infinity();
	for (unsigned i = 0; i < dimensions(); i++) {
		minimum = std::min(minimum, (*this)[i]);
	}
	return minimum;
}

double coordinate::euklidian(coordinate point) const{
	return coordinate(*this - point).length();
}

double coordinate::lengthSqr() const {
	double sumSqr = 0;
	for (unsigned i = 0; i < dimensions(); i++) {
		sumSqr += (*this)[i] * (*this)[i];
	}
	return sumSqr;
}

double coordinate::length() const {
	return sqrt(lengthSqr());
}

double coordinate::scp(const coordinate &a, const coordinate &b) {
	double scp = 0;
	for (unsigned i = 0; i < a.dimensions(); i++) {
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
