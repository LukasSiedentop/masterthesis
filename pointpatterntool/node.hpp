/*
 * node.hpp
 *
 *  Created on: 17.08.2015
 *      Author: lukas
 */

#ifndef NODE_HPP_
#define NODE_HPP_

#include <vector>

#include "coordinate.hpp"
#include "nodelist.hpp"

/**
 * Nodes of a pattern.
 */
class node {
private:
	// list containing this node. Only needed when pattern is periodic
	class nodelist* list;
	coordinate position;
	// TODO: shifter and neighbours as pair?! have to rewrite all loops over neighbours...
	// TODO: maybe store all possible shifters in the nodelist and only reference to them... saves memory and is more elegant
	std::vector<class node*> neighbours;
	std::vector<class coordinate> neighbourShifter;
	//std::vector<std::pair<class node*, class coordinate> > neighbours;
	// in the periodic case whether the node has links over the boundary, in the non-periodic case whether the node is farer away from the boundary than a given value
	bool edgenode;
public:
	node();
	node(class nodelist* list, coordinate pos);

	bool isEdgenode();
	coordinate getPosition() const;
	std::vector<class node*>* getNeighbours();
	std::vector<class coordinate>* getNeighbourShifters();
	bool isNeighbour(node* node);
	// sets the neighbourhood, together with a given shifter if the link lies over the boundary. Takes the shifters of this point and the given point and calculates from that the shifter of the link.
	void addNeighbour(node* n, coordinate ownShifter = coordinate(0, 0, 0), coordinate foreignShifter = coordinate(0, 0, 0));
	double euklidian(node* node);
	double euklidian(coordinate point);
	double euklidianPeriodic(node* node);
	double angle(node* nodeA, node* nodeB);
	double anglePeriodic(node* nodeA, node* nodeB);
	unsigned int countNeighbours();
	void shift(coordinate shifter);
	void scale(double a);
	void scaleAnisotropic(double ax, double ay, double az);

	// if this node is farther away than distance from the bounding box, it is not an edgenode
	void setEdgenode(double distance);

	// compares two nodes by comparing their position
	bool equals(node* node);

	// z-Metric
	static bool zmetric(const node* a, const node* b);
};

#endif /* NODE_HPP_ */
