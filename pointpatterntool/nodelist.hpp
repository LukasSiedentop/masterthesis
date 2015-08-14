/*
 * nodelist.hpp
 *
 * Datenstruktur des Punktmusters: eine doppelt verkettete Liste mit Kopf.
 *
 * TODO: Liste erbt von std::vector -> NodeHead = vector<Node>
 *
 * www.codeproject.com/Articles/668818/Implementing-a-Doubly-Linked-List-to-be-used-on-an
 * de.wikibooks.org/wiki/C%2B%2B-Programmierung
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

using namespace std;

/**
 * Kopf der Liste. Hier stecken Daten wie Listenlänge, erstes Element, letztes Element drin, die sonst jedesmal in jedem Knoten der Liste geändert werden müssten.
 */
class NodeHead {
private:
	/* Listenspezifisches */
	// Zeiger auf das erste und letzte Element
	class Node *first, *last;
	// Länge der Liste, mit 0 initialisiert
	int size;

	/* Datenspezifisches */
	double minX, minY, minZ, maxX, maxY, maxZ;
	// Periodizität des Musters
	bool periodic;
public:
	NodeHead();
	NodeHead(bool periodicity);
	~NodeHead();

	/* Listenspezifisches */
	// gibt die Listenlänge zurück
	int length();
	// fügt ein Element node am Ende der Liste an
	void add(Node * node);
	// löscht das Element node
	void remove(Node * node);
	// zeigt die Liste an.
	void display();
	// den ersten Knoten geben
	Node * getFirst();

	/* Datenspezifisches */
	// Minimalwerte setzen
	void setMins(double xMin, double yMin, double zMin);
	// Maximalwerte setzen
	void setMaxs(double xMax, double yMax, double zMax);
	// Bewegt die Liste um den gegebenen Vektor
	void shiftList(double sx, double sy, double sz);
	// Skaliert die Liste um die gegebenen Faktoren a,b, c
	void scaleList(double a, double b, double c);
	// Periodizität
	bool isPeriodic();
	// skaliert die Liste, sodass die Dichte der gewünschten entspricht
	void setDensity(double density);
	// berechnet die Punktdichte der Liste
	double getDensity();
	// überprüft ob ein Knoten an einem Punkt exisitiert
	Node * getAt(double x, double y, double z);
	// gibt die Minimalwerte der Box zurück
	double getMinX();
	double getMinY();
	double getMinZ();
	// gibt die Maximalwerte der Box zurück
	double getMaxX();
	double getMaxY();
	double getMaxZ();
	// gibt die Größe der Box zurück
	double lengthX();
	double lengthY();
	double lengthZ();
	// Gibt die 26 Vektoren zurück um das Muster periodisch fortzusetzen
	vector<vector<double> > getShifters();
	// Zählt die Punkte in einer gegebenen Kugel
	int pointsInside(double r, double mx, double my, double mz);
	// gibt Statistiken der Liste als String zurück TODO: was noch?
	string listStats(const char commentDelimeter[] = "\t");
};

/**
 * Element der Liste. Neben der Daten des Punktes im Punktmuster hat er drei Zeiger: das Element davor, dannach und den Kopf der Liste.
 */
class Node {
private:
	/* Listenspezifisches */
	// Zeiger auf den Listenkopf
	class NodeHead *head;
	// Zeiger auf Vorgänger und Nachfolger, mit NULL initialisiert
	Node *prev, *next;

	/* Datenspezifisches */
	// Position TODO: als array oder wie bei Dirk als extra Klasse
	double x, y, z;
	// Verbundene Nachbarn
	class Neighbour *neighbours;
	//Distanz zu Nachbarn, Winkel zu Nachbarn, Voronoi Nachbarn,...

	// prüft ob der gegebene Knoten schon ein Nachbar ist
	bool isNeighbour(Node * node);
public:
	Node();
	Node(double x, double y, double z);
	~Node();

	/* Listenspezifisches */
	// vertauscht diesen Knoten mit node.
	void swap(Node * node);
	// getter und setter
	void setHead(NodeHead * newHead);
	void setNext(Node * node);
	void setPrev(Node * node);
	void setX(double newX);
	NodeHead * getHead();
	Node * getNext();
	Node * getPrev();
	double getX();
	double getY();
	double getZ();
	Neighbour* getNeighbours();

	/* Datenspezifisches */
	// Gibt die Euklidsche Distanz zu node zurück.
	double euklidian(Node * node);
	// Gibt die Euklidsche Distanz zu einem gegebenen Punkt zurück.
	double euklidian(double px, double py, double pz);
	// Gibt die Euklidsche Distanz zu node zurück, mit periodische Randbedingungen angenommen
	double euklidianPeriodic(Node * node);
	// berechnet den Winkel zwischen diesem und den zwei gegebenen Knoten
	double angle(Node * nodeA, Node * nodeB);
	// berechnet den Winkel unter beachtung periodischer Randbedingungen
	double anglePeriodic(Node * nodeA, Node * nodeB);
	// fügt der Nachbarnliste einen weiteren Nachbarn hinzu
	void addNeighbour(Node * node);
	// vergleicht diesen Knoten mit dem gegebenen, wenn beide an der gleichen Position sind, sind sie gleich (TODO: Nachbarn und so auch vergleichen?)
	bool equals(Node * node);
	// zählt die Nachbarn
	int countNeighbours();
	// verschiebt den Knoten um den angegebenen Vektor
	void shift(double sx, double sy, double sz);
	// multipliziert die Knotenposition um die gegebenen Faktoren
	void scale(double a, double b, double c);
};

/**
 * Nachbarn eines Knotens in Form einer einfach verketteten Liste.
 */
class Neighbour {
private:
	Neighbour *next;
	Node * node;
public:
	Neighbour();
	Neighbour(Node * node);
	~Neighbour();

	// setter & getter
	Node * getNode();
	Neighbour * getNext();

	void add(Node * node);
};

#endif /* NODELIST_HPP_ */
