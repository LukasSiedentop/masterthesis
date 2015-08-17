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
 * Element der Liste. Neben der Daten des Punktes im Punktmuster hat er drei Zeiger: das Element davor, dannach und den Kopf der Liste.
 */
class node {
private:
	// Liste, zu der der Knoten gehört (gebraucht für periodizität)
	nodelist * list;
	// Position
	coordinate position;

	// Verbundene Nachbarn
	std::vector<class node *> neighbours;

	// Liste ob eine Koordinate über die grenze reicht (-1:in negativer Richtung drüber, 0: drin, +1: in positiver Richtung drüber)
	//std::vector<int> inBox;

	// prüft ob der gegebene Knoten schon ein Nachbar ist
	bool isNeighbour(node * node);
public:
	node();
	node(nodelist * list, double x, double y, double z);
	~node();

	coordinate getPosition();

	std::vector<class node *> * getNeighbours();
	//std::vector<int> * isInBox();

	/* Datenspezifisches */
	// Gibt die Euklidsche Distanz zu node zurück.
	double euklidian(node * node);
	// Gibt die Euklidsche Distanz zu einem gegebenen Punkt zurück.
	double euklidian(coordinate point);
	// Gibt die Euklidsche Distanz zu node zurück, mit periodische Randbedingungen angenommen
	double euklidianPeriodic(node * node);
	// berechnet den Winkel zwischen diesem und den zwei gegebenen Knoten
	double angle(node * nodeA, node * nodeB);
	// berechnet den Winkel unter beachtung periodischer Randbedingungen
	double anglePeriodic(node * nodeA, node * nodeB);
	// fügt der Nachbarnliste einen weiteren Nachbarn hinzu
	void addNeighbour(node * node);
	// vergleicht diesen Knoten mit dem gegebenen, wenn beide an der gleichen Position sind, sind sie gleich (TODO: Nachbarn und so auch vergleichen?)
	bool equals(node * node);
	// zählt die Nachbarn
	int countNeighbours();
	// verschiebt den Knoten um den angegebenen Vektor
	void shift(coordinate shifter);
	// multipliziert die Knotenposition um die gegebenen Faktoren
	void scale(double a, double b, double c);
};

#endif /* NODE_HPP_ */
