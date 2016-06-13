/*
 * node.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "node.hpp"

//using namespace std;

node::node() :
		position(coordinate()), edgenode(false) { //, inBox(NULL) {
	list = NULL;

	std::vector<class node*> neighs;
		neighbours = neighs;
		std::vector<class coordinate> neighShifters;
		neighbourShifter = neighShifters;
}

node::node(class nodelist* list, coordinate pos) :
		list(list), position(pos), edgenode(false) {	//, inBox(NULL) {
	std::vector<class node *> neighs;
	neighbours = neighs;

	std::vector<class coordinate> neighShifters;
	neighbourShifter = neighShifters;
}

void node::setEdgenode(double distance) {
	// TODO: this must be possible in less costly way!? hopefully without using the underlying list but its extend
	// Alle knoten die mehr als den 0-shifter haben sind randknoten
	edgenode =
			(list->getShifted(position - list->getMid(), distance).size() - 1);
}

bool node::isEdgenode() {
	return edgenode;
}

coordinate node::getPosition() const {
	return position;
}

std::vector<class coordinate>* node::getNeighbourShifters() {
	return &neighbourShifter;
}

std::vector<class node*>* node::getNeighbours() {
	return &neighbours;
}

double node::euklidian(node* node) {
	return euklidian(node->getPosition());
}

double node::euklidian(coordinate point) {
	return position.euklidian(point);
}

double node::euklidianPeriodic(node* node) {
	coordinate differenceVec = position - node->getPosition();

// Wenn die Länge größer als die Featuresize ist...
	double length = differenceVec.length();
	if (list->getMaxFeatureSize() < length) {
		differenceVec = coordinate::getVec(position, node->getPosition(),
				list->getShifters());
		length = differenceVec.length();
	}

	return length;
}

double node::angle(node * nodeA, node * nodeB) {
// cos alpha = skp(a,b) / |a|*|b|

	coordinate vec1 = nodeA->getPosition() - position;

	coordinate vec2 = nodeB->getPosition() - position;

	double len1sqr = vec1.lengthSqr();
	double len2sqr = vec2.lengthSqr();

	double skp = coordinate::scp(vec1, vec2);

	return acos(skp / sqrt(len1sqr * len2sqr));
}

double node::anglePeriodic(node* nodeA, node* nodeB) {
	coordinate vec1 = position - nodeA->getPosition();
	coordinate vec2 = position - nodeB->getPosition();

	double len1sqr = vec1.lengthSqr();
	double len2sqr = vec2.lengthSqr();

	double maxFeatureSize = list->getLengths().lengthSqr() / 4;

	if ((len1sqr > maxFeatureSize) || (len2sqr > maxFeatureSize)) {
		vec1 = coordinate::getVec(position, nodeA->getPosition(),
				list->getShifters());
		vec2 = coordinate::getVec(position, nodeB->getPosition(),
				list->getShifters());
		len1sqr = vec1.lengthSqr();
		len2sqr = vec2.lengthSqr();
	}

// Skalarprodukt
	double skp = coordinate::scp(vec1, vec2);

// Winkel ausgeben
	return acos(skp / sqrt(len1sqr * len2sqr));
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
	// TODO: valid?
	position *= coordinate(ax, ay, az);
	/*
	position[0] *= ax;
	position[1] *= ay;
	position[2] *= az;*/
}

void node::addNeighbour(node* n, coordinate ownShifter, coordinate foreignShifter ) {
	if (this->isNeighbour(n) || (this->equals(n))) {
		return;
	}
	neighbours.push_back(n);
	neighbourShifter.push_back(foreignShifter-ownShifter);
}

bool node::isNeighbour(node* node) {
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
