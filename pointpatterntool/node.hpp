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
	// list containing this node. Only needed if pattern is periodic
	class nodelist* list;
	coordinate position;
	std::vector<class node* > neighbours;
	bool edgenode;
public:
	node();
	node(class nodelist* list, coordinate pos);

	bool isEdgenode();
	coordinate getPosition() const;
	std::vector<class node* >* getNeighbours();
	bool isNeighbour(node* node);
	void addNeighbour(node* n);
	double euklidian(node* node);
	double euklidian(coordinate point);
	double euklidianPeriodic(node* node);
	double angle(node* nodeA, node* nodeB);
	double anglePeriodic(node* nodeA, node* nodeB);
	int countNeighbours();
	void shift(coordinate shifter);
	void scale(double a);

	// if this node is farther away than distance from the bounding box, it is not an edgenode
	void setEdgenode(double distance);

	// compares two nodes by comparing their position
	bool equals(node* node);
};

#endif /* NODE_HPP_ */
