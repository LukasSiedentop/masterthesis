/*
 * nodelist.hpp
 *
 * Datenstruktur des Punktmusters: eine doppelt verkettete Liste mit Kopf.
 *
 * www.codeproject.com/Articles/668818/Implementing-a-Doubly-Linked-List-to-be-used-on-an
 * de.wikibooks.org/wiki/C%2B%2B-Programmierung
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

#ifndef NODELIST_HPP_
#define NODELIST_HPP_

#include <iostream>
#include <cmath>


/**
 * Kopf der Liste. Hier stecken Daten wie Listenlänge, erstes Element, letztes Element drin, die sonst jedesmal in jedem Knoten der Liste geändert werden müssten.
 */
class NodeHead {
private:
	/* Listenspezifisches */
	class Node * first, *last;
	int size;

	/* Datenspezifisches */
	//... ?
public:
	NodeHead();
	~NodeHead();

	/* Listenspezifisches */
	// gibt die Listenlänge zurück
	int length();
	// fügt ein Element node am Ende der Liste an
	void add(Node * node);
	// löscht das Element node
	void remove(Node * node);
	//...

	/* Datenspezifisches */
	//... ?
};

/**
 * Element der Liste. Neben der Daten des Punktes im Punktmuster hat er drei Zeiger: das Element davor, dannach und den Kopf der Liste.
 */
class Node {
private:
	/* Listenspezifisches */
	class NodeHead * head;
	Node * prev, *next;


	/* Datenspezifisches */
	// Position
	double x, y, z;
	// Verbundene Nachbarn, Distanz zu Nachbarn, Winkel zu Nachbarn, Voronoi Nachbarn,...
public:
	Node();
	Node(double x, double y, double z);
	~Node();

	/* Listenspezifisches */
	// vertauscht diesen Knoten mit node.
	void swap(Node * node);
	// getter und setter (warum nochmal? brauchen wir wirklich private und so?)
	void setHead(NodeHead * newHead);
	void setNext(Node * node);
	void setPrev(Node * node);
	void setX(double newX);
	NodeHead * getHead();
	Node * getNext();
	Node * getPrev();
	double getX();
	//...

	/* Datenspezifisches */
	// Gibt die Euklidsche Distanz zu node zurück.
	double euklidian(Node * node);
};


#endif /* NODELIST_HPP_ */
