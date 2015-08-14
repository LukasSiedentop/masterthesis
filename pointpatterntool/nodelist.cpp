/*
 * nodelist.cpp
 *
 *  Created on: 24.07.2015
 *      Author: Lukas Siedentop
 */

#include "nodelist.hpp"

//#include <iostream>
#include <math.h>
#include <sstream>


using namespace std;

/* NodeHead */

NodeHead::NodeHead() :
		first(NULL), last(NULL), size(0), minX(0), minY(0), minZ(0), maxX(0), maxY(
				0), maxZ(0), periodic(0) {
}

NodeHead::NodeHead(bool periodicity) :
		first(NULL), last(NULL), size(0), minX(0), minY(0), minZ(0), maxX(0), maxY(
				0), maxZ(0), periodic(periodicity) {
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
		first = last = node;
	}

}

Node * NodeHead::getFirst() {
	return first;
}

void NodeHead::setMins(double xMin, double yMin, double zMin) {
	minX = xMin;
	minY = yMin;
	minZ = zMin;
}

void NodeHead::setMaxs(double xMax, double yMax, double zMax) {
	maxX = xMax;
	maxY = yMax;
	maxZ = zMax;
}

void NodeHead::shiftList(double sx, double sy, double sz) {
	minX += sx;
	minY += sy;
	minZ += sz;

	maxX += sx;
	maxY += sy;
	maxZ += sz;

	Node * nodeIter = first;
	while (nodeIter) {
		nodeIter->shift(sx, sy, sz);
		nodeIter = nodeIter->getNext();
	}
}

void NodeHead::scaleList(double a, double b, double c) {
	Node * nodeIter = first;
	while (nodeIter) {
		nodeIter->scale(a, b, c);
		nodeIter = nodeIter->getNext();
	}
}

bool NodeHead::isPeriodic() {
	return periodic;
}

void NodeHead::setDensity(double density) {
	double factor = density / getDensity();
	scaleList(factor, factor, factor);
}

double NodeHead::getDensity() {
	return size / (lengthX() * lengthY() * lengthZ());
}

Node * NodeHead::getAt(double x, double y, double z) {

	// Toleranz
	double t = 0.000001;

	// Knoteniteration
	Node* node = first;
	while (node) {

		double xDist = fabs(x - node->getX());
		double yDist = fabs(y - node->getY());
		double zDist = fabs(z - node->getZ());

		// Minimalen Abstand finden, wenn das Muster periodisch fortsetzbar ist
		if (periodic) {
			xDist = min(xDist, fabs(xDist - lengthX()));
			yDist = min(yDist, fabs(yDist - lengthY()));
			zDist = min(zDist, fabs(zDist - lengthZ()));
		}

		// Distanzen mit der Toleranz vergleichen
		if ((xDist < t && yDist < t && zDist < t)) {
			return node;
		}
		node = node->getNext();
	}

	return NULL;
}

double NodeHead::getMinX() {
	return minX;
}

double NodeHead::getMinY() {
	return minY;
}

double NodeHead::getMinZ() {
	return minZ;
}

double NodeHead::getMaxX() {
	return maxX;
}

double NodeHead::getMaxY() {
	return maxY;
}

double NodeHead::getMaxZ() {
	return maxZ;
}

double NodeHead::lengthX() {
	return maxX - minX;
}

double NodeHead::lengthY() {
	return maxY - minY;
}

double NodeHead::lengthZ() {
	return maxZ - minZ;
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
			neighboursIter = neighboursIter->getNext();
		}
		std::cout << std::endl;

		nodesIter = nodesIter->getNext();
		i++;
	}

}

vector<vector<double> > NodeHead::getShifters() {
	vector<vector<double> > shifters;

	// Boxlängen
	double lx = lengthX(), ly = lengthY(), lz = lengthZ();

	// alle Kombinationen durchgehen
	for (int ix = -1; ix < 2; ix++) {
		for (int iy = -1; iy < 2; iy++) {
			for (int iz = -1; iz < 2; iz++) {
				vector<double> coord;
				coord.push_back(ix * lx);
				coord.push_back(iy * ly);
				coord.push_back(iz * lz);

				shifters.push_back(coord);
			}
		}
	}

	return shifters;
}

int NodeHead::pointsInside(double r, double mx, double my, double mz) {
	Node * nodeIter = first;

	int ctr = 0;
	vector<vector<double> > shifters = getShifters();
	// TODO: paralelisieren
	while (nodeIter) {

		if (periodic) {
			// Abstandsvergleich mit in jeder Raumrichtung verschobenen Kugel
			for (unsigned int i = 0; i < shifters.size(); i++) {
				if (r
						> nodeIter->euklidian(mx + shifters[i][0],
								my + shifters[i][1], mz + shifters[i][2])) {
					ctr++;
				}
			}
		} else if (r > nodeIter->euklidian(mx, my, mz)) {
			ctr++;
		}

		nodeIter = nodeIter->getNext();
	}

	return ctr;
}

string NodeHead::listStats(const char commentDelimeter[]) {
	double volume = lengthX()*lengthY()*lengthZ();
	double density = getDensity();

	stringstream stream;
	stream << commentDelimeter << "Anzahl Knoten: " << size << endl;
	stream << commentDelimeter << "Extremalwerte: (" << minX << ", " << minY << ", " << minZ << "), (" << maxX << ", " << maxY << ", " << maxZ << ")" << endl;
	stream << commentDelimeter << "Volumen: " << volume << endl;
	stream << commentDelimeter << "Punktdichte: " << density << endl;

	return stream.str();
}

/* Node */

Node::Node() :
		head(NULL), prev(NULL), next(NULL), x(0), y(0), z(0), neighbours(NULL) {
}

Node::Node(double x, double y, double z) :
		head(NULL), prev(NULL), next(NULL), x(x), y(y), z(z), neighbours(NULL) {
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

void Node::setX(double newX) {
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
	return euklidian(node->getX(), node->getY(), node->getZ());
}

double Node::euklidian(double px, double py, double pz) {

	double xDist = x - px;
	double yDist = y - py;
	double zDist = z - pz;

	return sqrt(pow(xDist, 2) + pow(yDist, 2) + pow(zDist, 2));
}

double Node::euklidianPeriodic(Node * node) {
	// Minimalen Abstand finden, wenn das Muster periodisch fortsetzbar ist
	double xDist = min(fabs(x - node->getX()),
			fabs(fabs(x - node->getX()) - head->lengthX()));
	double yDist = min(fabs(y - node->getY()),
			fabs(fabs(y - node->getY()) - head->lengthY()));
	double zDist = min(fabs(z - node->getZ()),
			fabs(fabs(z - node->getZ()) - head->lengthZ()));

	return sqrt(pow(xDist, 2) + pow(yDist, 2) + pow(zDist, 2));
}

double Node::angle(Node * nodeA, Node * nodeB) {
	// cos alpha = skp(a,b) / |a|*|b|

	double vec1X = nodeA->getX() - x;
	double vec1Y = nodeA->getY() - y;
	double vec1Z = nodeA->getZ() - z;

	double vec2X = nodeB->getX() - x;
	double vec2Y = nodeB->getY() - y;
	double vec2Z = nodeB->getZ() - z;

	double len1 = pow(vec1X, 2) + pow(vec1Y, 2) + pow(vec1Z, 2);
	double len2 = pow(vec2X, 2) + pow(vec2Y, 2) + pow(vec2Z, 2);

	double skp = vec1X * vec2X + vec1Y * vec2Y + vec1Z * vec2Z;

	return acos(skp / sqrt(len1 * len2));
}

double Node::anglePeriodic(Node * nodeA, Node * nodeB) {
	double vec1X = x - nodeA->getX();
	double vec1Y = y - nodeA->getY();
	double vec1Z = z - nodeA->getZ();

	double vec2X = x - nodeB->getX();
	double vec2Y = y - nodeB->getY();
	double vec2Z = z - nodeB->getZ();

	double len1 = pow(vec1X, 2) + pow(vec1Y, 2) + pow(vec1Z, 2);
	double len2 = pow(vec2X, 2) + pow(vec2Y, 2) + pow(vec2Z, 2);

	// Wenn Vektor 1 länger als die halbe Boxlänge ist...
	if (sqrt(len1) > head->lengthX() / 2) {
		// temporär die Distanz zwischenspeichern
		double dist = 0;
		vector<vector<double> > shifters = head->getShifters();
		for (unsigned int i = 0; i < shifters.size(); i++) {
			// Distanz berechnen
			// wenn deren betrag kleiner ist als die vektorkomponente, austauschen

			// Vektor A
			dist = x - (nodeA->getX() + shifters[i][0]);
			if (fabs(dist) < fabs(vec1X)) {
				vec1X = dist;
			}

			dist = y - (nodeA->getY() + shifters[i][1]);
			if (fabs(dist) < fabs(vec1Y)) {
				vec1Y = dist;
			}

			dist = z - (nodeA->getZ() + shifters[i][2]);
			if (fabs(dist) < fabs(vec1Z)) {
				vec1Z = dist;
			}
		}
		// Länge neu berechnen
		len1 = pow(vec1X, 2) + pow(vec1Y, 2) + pow(vec1Z, 2);
	}

	// Wenn Vektor 2 länger als die halbe Boxlänge ist...
	if (sqrt(len2) > head->lengthX() / 2) {
		// temporär die Distanz zwischenspeichern
		double dist = 0;
		vector<vector<double> > shifters = head->getShifters();
		for (unsigned int i = 0; i < shifters.size(); i++) {
			// Distanz berechnen
			// wenn deren betrag kleiner ist als die vektorkomponente, austauschen

			// Vektor B
			dist = x - (nodeB->getX() + shifters[i][0]);
			if (fabs(dist) < fabs(vec2X)) {
				vec2X = dist;
			}

			dist = y - (nodeB->getY() + shifters[i][1]);
			if (fabs(dist) < fabs(vec2Y)) {
				vec2Y = dist;
			}

			dist = z - (nodeB->getZ() + shifters[i][2]);
			if (fabs(dist) < fabs(vec2Z)) {
				vec2Z = dist;
			}
		}
		// Länge neu berechnen
		len2 = pow(vec2X, 2) + pow(vec2Y, 2) + pow(vec2Z, 2);
	}

	// Skalarprodukt
	double skp = vec1X * vec2X + vec1Y * vec2Y + vec1Z * vec2Z;

	// Winkel ausgeben
	return acos(skp / sqrt(len1 * len2));
}

int Node::countNeighbours() {
	int num = 0;

	Neighbour * neighIter = neighbours;

	while (neighIter) {
		num++;
		neighIter = neighIter->getNext();
	}

	return num;
}

void Node::shift(double sx, double sy, double sz) {
	x += sx;
	y += sy;
	z += sz;
}

void Node::scale(double a, double b, double c) {
	x *= a;
	y *= b;
	z *= c;
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
		neighIter = neighIter->getNext();
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

/* Neighbour */

Neighbour::Neighbour() :
		next(NULL), node(NULL) {
}

Neighbour::Neighbour(Node * node) :
		next(NULL), node(node) {
}

Neighbour * Neighbour::getNext() {
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
