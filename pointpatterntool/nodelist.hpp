/*
 * nodelist.hpp
 *
 *  Created on: 24.07.2015
 *      Author: Lukas Siedentop
 */

#ifndef NODELIST_HPP_
#define NODELIST_HPP_

#include <iostream>
#include <cmath>
#include <math.h>
#include <sstream>
#include <vector>
#include <boost/progress.hpp>
#include <ctime>
#include <iterator>
#include <memory>

#include "coordinate.hpp"
#include "node.hpp"
#include "functions.hpp"

using namespace std;

/**
 * Datastructure to hold points of a pattern. Methods for analysis (like distance to neighbours distribution, hyperuniformity,...) of the pattern are provided here.
 */
class nodelist {
private:
	bool periodic;
	string name;
	// Thou shalt not inherit from STL
	std::vector<class node*> list;
	// Bounding box
	coordinate min, max;

	coordinate getMins();
	coordinate getMaxs();
	bool isPeriodic() const;
	void setDensity(double density);
	double getDensity();
	double getVolume();
	int countEdgenodes();

	// sets the neighbours of each point to the closest four (or less) points
	void setNeighbours();
	// counts the point within a given sphere
	int pointsInside(const vector<coordinate>& points, const coordinate& mid, const double r) const;

public:
	nodelist();
	// constructs empty list
	nodelist(bool periodicity, string name);
	// constructs a pattern (density of points 1, within 10^3 cubicle) with: pattern=1 - random points, pattern=2 - points arranged in a diamond lattice.
	nodelist(int pattern, bool periodicity);

	std::vector<node*>::iterator begin();
	std::vector<node*>::iterator end();
	void setEdgenodes(double distance);
	coordinate getLengths();
	coordinate getMid();
	double getMaxFeatureSize();
	void setMins(coordinate mins);
	void setMaxs(coordinate maxs);
	vector<double> neighbourDistribution();
	vector<double> lengthDistribution();
	vector<double> angleDistribution();
	vector<vector<double> > hyperuniformity(unsigned int nr=50, unsigned int n=100);
	string getName();
	void shiftList(coordinate shifter);
	void scaleList(double a);
	void deleteEntries();

	// subscript operator
	//const node& operator[](const int i) const;
	node& operator[](int i);

	// returns the number of nodes
	int size();
	// Adds a node at the given coordinate if none exists there and returns it in order to be able to set the neighbourhood.
	node* add(double x, double y, double z);
	// Returns the pattern in a matrix gnuplot can interpret.
	vector<vector<double> > getGnuplotMatrix();
	// Returns the 26 shifting vectors to continue the pattern periodically. TODO: n-Dimensional
	vector<coordinate> getShifters();
	// Returns the 26 shifting vectors to continue the pattern periodically. Assumes the pattern is around (0,0,0). Only necessary vectors are returned so that the given box lies within the continued pattern. TODO: n-Dimensional
	vector<coordinate> getShifted(coordinate mid, double halfExtend);
	// Writes the pattern in a representation of a povray scene. The radius r of the cylinders if determined while rendering.
	void writePOV();
	// normalizes the pattern to density of points=1, midpoint of bounding box = (0,0,0)
	double normalize();
	// returns statistics of the list
	string listStats(const char commentDelimeter[] = "\t");
	// outputs list to console
	void display();

};

#endif /* NODELIST_HPP_ */
