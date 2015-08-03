/*
 * nodelist.cpp
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

<<<<<<< HEAD
#include "nodelist.hpp"

NodeHead::NodeHead() :
		first(NULL), last(NULL), size(-1) {
=======
#include <math.h>
#include "nodelist.hpp"

NodeHead::NodeHead() :
		first(NULL), last(NULL), size(0), minX(0), minY(0), minZ(0), maxX(0), maxY(0), maxZ(0) {
>>>>>>> 8c44c62822241421bedf5ec48c67ff453b159266
}

int NodeHead::length() {
	return size;
}

void NodeHead::add(Node * node) {
	// Kopf setzen
	node->setHead(this);

	// Größe anpassen
	size++;

	if (first) { // Liste ist nicht leer
		// Zeiger umbiegen
		last->setNext(node);
		node->setPrev(last);
		last = node;
	} else { // Liste ist leer
		// Zeiger setzen
<<<<<<< HEAD
		first = node;
		last = node;
	}

=======
		first = last = node;
	}

}

Node * NodeHead::getFirst() {
	return first;
}

Node * NodeHead::getAt(double x, double y, double z) {

	// Toleranz
	double t = 0.000001;

	Node* node = first;

	while (node) {
		// TODO: reicht der direkte Vergleich immer? Vielleicht einen kleinen Radius definieren, wie Dirk
		if ((fabs(x - node->getX()) < t && fabs(y - node->getY()) < t
				&& fabs(z - node->getZ()) < t)) {
			return node;
		}
		node = node->getNext();
	}

	return NULL;
}

void NodeHead::display() {
	std::cout << "Anzahl Elemente: " << size << std::endl;
	Node * nodesIter = first;
	int i = 1;
	while (nodesIter) {
		std::cout << char(9) << i << "-tes Element: (" << nodesIter->getX()
				<< ", " << nodesIter->getY() << ", " << nodesIter->getZ()
				<< ") mit " << nodesIter->countNeighbours()
				<< " Nachbarn bei: ";
		Neighbour* neighboursIter = nodesIter->getNeighbours();
		while (neighboursIter) {
			std::cout << "(" << neighboursIter->getNode()->getX() << ", "
					<< neighboursIter->getNode()->getY() << ", "
					<< neighboursIter->getNode()->getZ() << ") ";
			neighboursIter = neighboursIter->getNextNeighbour();
		}
		std::cout << std::endl;

		nodesIter = nodesIter->getNext();
		i++;
	}
>>>>>>> 8c44c62822241421bedf5ec48c67ff453b159266

}

Node::Node() :
<<<<<<< HEAD
		head(NULL), prev(NULL), next(NULL), x(0), y(0), z(0) {
}

Node::Node(double x, double y, double z) :
		head (NULL), prev(NULL), next(NULL), x(x), y(y), z(z) {
=======
		head(NULL), prev(NULL), next(NULL), x(0), y(0), z(0), neighbours(NULL) {
}

Node::Node(double x, double y, double z) :
		head(NULL), prev(NULL), next(NULL), x(x), y(y), z(z), neighbours(NULL) {
>>>>>>> 8c44c62822241421bedf5ec48c67ff453b159266
}

void Node::setHead(NodeHead * newHead) {
	head = newHead;
}

void Node::setNext(Node * node) {
	next = node;
}

void Node::setPrev(Node * node) {
	prev = node;
}

<<<<<<< HEAD
void Node::setX(double newX){
=======
void Node::setX(double newX) {
>>>>>>> 8c44c62822241421bedf5ec48c67ff453b159266
	x = newX;
}

NodeHead * Node::getHead() {
	return head;
}

Node * Node::getNext() {
	return next;
}

Node * Node::getPrev() {
	return prev;
}

<<<<<<< HEAD
double Node::getX(){
	return x;
}
=======
double Node::getX() {
	return x;
}

double Node::getY() {
	return y;
}

double Node::getZ() {
	return z;
}

Neighbour * Node::getNeighbours() {
	return neighbours;
}

double Node::euklidian(Node * node) {
	return sqrt(
			pow((x - node->getX()), 2) + pow((y - node->getY()), 2)
					+ pow((z - node->getZ()), 2));
}

double Node::angle(Node * nodeA, Node * nodeB) {
	// cos alpha = |a*b| / |a|*|b|

	double vec1X = x-nodeA->getX();
	double vec1Y = y-nodeA->getY();
	double vec1Z = z-nodeA->getZ();
	double len1 = pow(vec1X,2) + pow(vec1Y,2) + pow(vec1Z,2);

	double vec2X = x-nodeB->getX();
	double vec2Y = y-nodeB->getY();
	double vec2Z = z-nodeB->getZ();
	double len2 = pow(vec2X,2) + pow(vec2Y,2) + pow(vec2Z,2);

	double skp = vec1X*vec2X + vec1Y*vec2Y + vec1Z*vec2Z;

	return acos(sqrt(skp/(len1*len2)));
}

int Node::countNeighbours() {
	int num = 0;

	Neighbour * neighIter = neighbours;

	while (neighIter) {
		num++;
		neighIter = neighIter->getNextNeighbour();
	}

	return num;
}

void Node::addNeighbour(Node * node) {
	// aufhören wenn node schon ein Nachbar ist bzw node dieser Knoten ist
	if (this->isNeighbour(node) || (this->equals(node))) {
		return;
	}

	// wenn es schon einen Nachbarn zu diesem Knoten gibt...
	if (neighbours) {
		// ...diesen den neuen hinzufügen lassen
		neighbours->add(node);
	} else {
		// ...ansonsten initialisieren.
		neighbours = new Neighbour(node);
	}
}

bool Node::isNeighbour(Node * node) {
	// iteriert durch alle Nachbarn und vergleicht die Positionen der Knoten
	Neighbour *neighIter = neighbours;
	if (neighIter) {
		// TODO: besser ist ein wirklicher Vergleich der Speicherpositionen, nachschauen. sowas wie Pointer1 == Pointer2
		if (node->equals(neighIter->getNode())) {
			// ein gleicher ist gefunden
			return true;
		}
		neighIter = neighIter->getNextNeighbour();
	}

	// wenn kein gleicher gefunden ist, ists kein Nachbar
	return false;
}

bool Node::equals(Node * node) {
	// node == this, wenn die Positionen gleich sind.
	// Toleranz
	double t = 0.000001;
	return (fabs(x - node->getX()) < t && fabs(y - node->getY()) < t
			&& fabs(z - node->getZ()) < t);
}

Neighbour::Neighbour() :
		next(NULL), node(NULL) {
}

Neighbour::Neighbour(Node * node) :
		next(NULL), node(node) {
}

Neighbour * Neighbour::getNextNeighbour() {
	return next;
}

Node * Neighbour::getNode() {
	return node;
}

void Neighbour::add(Node * newNode) {
	// wenn der neue Knoten mit diesem übereinstimmt, abbrechen
	if (newNode->equals(node)) {
		return;
	}

	// wenn dieser Nachbar einen nächsten Nachbarn hat....
	if (next) {
		// ...den Nächsten Nachbarn den Knoten hinzufügen lassen
		next->add(newNode);
	} else {
		// ...ansonsten node als Nachbarn hinzufügen
		next = new Neighbour(newNode);
	}
}
>>>>>>> 8c44c62822241421bedf5ec48c67ff453b159266
