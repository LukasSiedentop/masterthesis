/*
 * node.cpp
 *
 *  Created on: 14.08.2015
 *      Author: Lukas Siedentop
 */

#include "node.hpp"

node::node() :
		position(coordinate()), edgenode(false) {
	list = NULL;

	std::vector<class node*> neighs;
		neighbours = neighs;
		std::vector<class coordinate> neighShifters;
		neighbourShifter = neighShifters;
}

node::node(class nodelist* list, coordinate pos) :
		list(list), position(pos), edgenode(false) {
	std::vector<class node *> neighs;
	neighbours = neighs;

	std::vector<class coordinate> neighShifters;
	neighbourShifter = neighShifters;
}

void node::setEdgenode(double distance) {
	edgenode = !position.insideAABB(list->getMins()+distance, list->getMaxs()-distance);
}

bool node::isEdgenode() const {
	return edgenode;
}

coordinate node::getPosition() const {
	return position;
}

std::vector<class coordinate>* node::getNeighbourShifters() {
	return &neighbourShifter;
}

coordinate* node::getNeighbourShifter(const unsigned int i) {
	return &neighbourShifter[i];
}

coordinate node::getShiftedNeighbourposition(const unsigned int& i) const{
	return neighbours[i]->getPosition() - neighbourShifter[i];
}

std::vector<class node*>* node::getNeighbours() {
	return &neighbours;
}

node* node::getNeighbour(const unsigned int i) {
	return neighbours[i];
}

unsigned int node::countNeighbours() {
	if (neighbours.size()!=neighbourShifter.size()){
		std::cout << "Error: not as many shifters as neighbours!" << std::endl;
	}
	return neighbours.size();
}

void node::shift(coordinate shifter) {
	position += shifter;
}

void node::scale(double a) {
	position *= a;
}

void node::scaleAnisotropic(double ax, double ay, double az) {
	position *= coordinate(ax, ay, az);
}

void node::addNeighbour(node* n, coordinate ownShifter, coordinate foreignShifter ) {
	// TODO: is not yet fully elegant: the neighbours are added uniquely, but the fact that one neighbour can have multiple shifters is neglected. Should only be relevant for small point densities.

	if (this->isNeighbour(n) || (this->equals(n))) {
		return;
	}
	neighbours.push_back(n);
	// subtract own shifter, as it is done so when adding the node to the list
	neighbourShifter.push_back(foreignShifter-ownShifter);
}

bool node::isNeighbour(node* node){
	for (std::vector<class node*>::iterator it = neighbours.begin();
			it != neighbours.end(); ++it) {
		if ((*it) == node) {
			return true;
		}
	}
	return false;
}

bool node::equals(node* node) {
	// compare positions
	return (position == node->getPosition());
}

bool node::zmetric(const node* a, const node* b) {
	// TODO: generalizable with nmetric, n being the dimension to compare
	return a->getPosition().z() < b->getPosition().z();
}
