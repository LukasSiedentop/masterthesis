/*
 * coordinate.cpp
 *
 *  Created on: 14.08.2015
 *      Author: Lukas Siedentop
 */

#include "coordinate.hpp"

coordinate::coordinate() {
}

coordinate::coordinate(unsigned int n) {
	position.resize(n);
	for (unsigned int i = 0; i < n; i++) {
		position[i] = (double) rand() / RAND_MAX;
	}
}

coordinate::coordinate(double pos, unsigned int n) {
	position.resize(n);
	for (unsigned int i = 0; i < n; i++) {
		position[i] = (pos);
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

std::vector<double>* coordinate::getVector(){
	return &position;
}

unsigned int coordinate::dimensions() const {
	return position.size();
}

double& coordinate::operator [](const int i) {
	return position[i];
}

const double& coordinate::operator [](const int i) const {
	return position[i];
}

bool coordinate::operator ==(const coordinate &rhs) const {
	for (unsigned int i = 0; i < dimensions(); i++) {
		if (std::abs(rhs[i] - position[i]) > tolerance) {
			return false;
		}
	}
	return true;
}

bool coordinate::operator !=(const coordinate &rhs) const {
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
		(*this)[i] += summand;
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

std::ostream& operator <<(std::ostream &os, const coordinate &obj) {
	return os << obj.toString();
}

std::string coordinate::toString(const std::string begin,
		const std::string delimiter, const std::string end,
		const int precision) const {
	std::stringstream stream;

	// write precision=4 non-zero digits for each double
	stream.precision(precision);

	stream << begin;

	for (unsigned i = 0; i < dimensions() - 1; i++) {
		stream << position[i] << delimiter;
	}

	stream << position[dimensions() - 1] << end;

	return stream.str();
}

bool coordinate::insideAABB(const coordinate& min, const coordinate& max) const {
	if (x()>max.x()) return false;
	if (x()<min.x()) return false;
	if (y()>max.y()) return false;
	if (y()<min.y()) return false;
	if (z()>max.z()) return false;
	if (z()<min.z()) return false;
	return true;
}
bool coordinate::insideAABB(const coordinate& mid, const double& r) const {
	return (fabs(x() - mid[0]) < r) && (fabs(y() - mid[1]) < r)
			&& (fabs(z() - mid[2]) < r);
}

double coordinate::x() const {
	return position[0];
}
double coordinate::y() const {
	return position[1];
}
double coordinate::z() const {
	return position[2];
}

double coordinate::min() const {
	double minimum = std::numeric_limits<double>::infinity();
	for (unsigned i = 0; i < dimensions(); i++) {
		minimum = std::min(minimum, position[i]);
	}
	return minimum;
}

double coordinate::max() const {
	double maximum = -std::numeric_limits<double>::infinity();
	for (unsigned i = 0; i < dimensions(); i++) {
		maximum = std::max(maximum, position[i]);
	}
	return maximum;
}

double coordinate::euklidian(const coordinate& point) const {
	return coordinate(*this - point).length();
}

double coordinate::lengthSqr() const {
	double sumSqr = 0;
	for (unsigned i = 0; i < dimensions(); i++) {
		sumSqr +=position[i] * position[i];
	}
	return sumSqr;
}

double coordinate::length() const {
	return sqrt(lengthSqr());
}

double coordinate::angle(const coordinate& b) const {
	return acos(this->scp(b) / sqrt(this->lengthSqr() * b.lengthSqr()));
}

std::vector<coordinate> coordinate::getShifters() const {
	std::vector<coordinate> shifters;

	// go through all combinations
	for (int ix = -1; ix < 2; ix++) {
		for (int iy = -1; iy < 2; iy++) {
			for (int iz = -1; iz < 2; iz++) {
				// add the shifting vector
				shifters.push_back(coordinate(ix * x(), iy * y(), iz * z()));
			}
		}
	}
	return shifters;
}

double coordinate::scp(const coordinate &b) const{
	double scp = 0;
	for (unsigned i = 0; i < this->dimensions(); i++) {
		scp += position[i] * b[i];
	}
	return scp;
}

coordinate coordinate::cpr(const coordinate &b) const{
	return coordinate((position[1] * b[2] - position[2] * b[1]), (position[2] * b[0] - position[0] * b[2]),
			(position[0] * b[1] - position[1] * b[0]));
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
