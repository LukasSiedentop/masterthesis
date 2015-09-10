/*
 * node.cpp
 *
 *  Created on: 14.08.2015
 *      Author: lukas
 */

#include "node.hpp"

using namespace std;

node::node() :
		position(coordinate()), edgenode(false) { //, inBox(NULL) {
	// TODO: geht das besser?
	nodelist* l = new nodelist();
	list = l;

	std::vector<class node*> neighs;
	neighbours = neighs;
}

node::node(class nodelist* list, coordinate pos) :
		list(list), position(pos), edgenode(false) {	//, inBox(NULL) {
	std::vector<class node *> neighs;
	neighbours = neighs;
}

void node::setList(nodelist* l) {
	list = l;
}

void node::setEdgenode(double distance) {
	// Alle knoten die mehr als den 0-shifter haben sind randknoten
	edgenode =
			(list->getShifted(position - list->getMid(), distance).size() - 1);
}

bool node::isEdgenode() {
	return edgenode;
}

coordinate* node::getPosition() const {
	return new coordinate(position);
}

std::vector<class node*>* node::getNeighbours() {
	return &neighbours;
}

double node::euklidian(node* node) {
	return euklidian(*node->getPosition());
}

double node::euklidian(coordinate point) {
	return position.euklidian(point);
}

double node::euklidianPeriodic(node* node) {
	coordinate differenceVec = position - *node->getPosition();

// Wenn die Länge größer als die Featuresize ist...
	double length = differenceVec.length();
	if (list->getMaxFeatureSize() < length) {
		differenceVec = coordinate::getVec(position, *node->getPosition(),
				list->getShifters());
		length = differenceVec.length();
	}

	return length;
}

double node::angle(node * nodeA, node * nodeB) {
// cos alpha = skp(a,b) / |a|*|b|

	coordinate vec1 = *nodeA->getPosition() - position;

	coordinate vec2 = *nodeB->getPosition() - position;

	double len1sqr = vec1.lengthSqr();
	double len2sqr = vec2.lengthSqr();

	double skp = coordinate::scp(vec1, vec2);

	return acos(skp / sqrt(len1sqr * len2sqr));
}

double node::anglePeriodic(node* nodeA, node* nodeB) {
	coordinate vec1 = position - *nodeA->getPosition();
	coordinate vec2 = position - *nodeB->getPosition();

	double len1sqr = vec1.lengthSqr();
	double len2sqr = vec2.lengthSqr();

	double maxFeatureSize = list->getLengths().lengthSqr() / 4;

	if ((len1sqr > maxFeatureSize) || (len2sqr > maxFeatureSize)) {
		vec1 = coordinate::getVec(position, *nodeA->getPosition(),
				list->getShifters());
		vec2 = coordinate::getVec(position, *nodeB->getPosition(),
				list->getShifters());
		len1sqr = vec1.lengthSqr();
		len2sqr = vec2.lengthSqr();
	}

// Skalarprodukt
	double skp = coordinate::scp(vec1, vec2);

// Winkel ausgeben
	return acos(skp / sqrt(len1sqr * len2sqr));
}

int node::countNeighbours() {
	return neighbours.size();
}

void node::shift(coordinate shifter) {
	position += shifter;
}

void node::scale(double a) {
	position *= a;
}

void node::addNeighbour(node* node) {
// aufhören wenn node schon ein Nachbar ist bzw node dieser Knoten ist
	if (this->isNeighbour(node) || (this->equals(node))) {
		return;
	}

	neighbours.push_back(node);
}

bool node::isNeighbour(node* node) {

	bool exists = 0;
	for (vector<class node*>::iterator it = neighbours.begin();
			it != neighbours.end(); ++it) {

		// TODO: geht das wirklich?
		// TODO: abbrechen wenn er existiert
		exists += ((*it) == node);
	}
	return exists;
}

bool node::equals(node* node) {
	/*
	 // Nachbarn vergleichen
	 bool equals = 0;
	 for (vector<node *>::iterator neighbour = neighbours.begin();
	 neighbour != neighbours.end(); ++neighbour) {
	 for (vector<node *>::iterator nodeNeighbour = neighbours.begin();
	 nodeNeighbour != neighbours.end(); ++nodeNeighbour) {

	 // TODO: prüfen
	 equals += (neighbour == nodeNeighbour);

	 }
	 }

	 //return ((position == node->getPosition()) && equals);
	 */

// Positionen vergleichen
	return (position == *node->getPosition());
}

bool node::operator <=(const node &rhs) {
	return position.lengthSqr() <= rhs.getPosition()->lengthSqr();
}
