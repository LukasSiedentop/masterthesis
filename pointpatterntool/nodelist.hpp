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

#include <boost/progress.hpp>

#include <ctime>

#include "coordinate.hpp"
#include "node.hpp"
#include "functions.hpp"

using namespace std;

/**
 * Kopf der Liste. Hier stecken Daten wie Listenlänge, erstes Element, letztes Element drin, die sonst jedesmal in jedem Knoten der Liste geändert werden müssten.
 */
class nodelist: public std::vector<class node*> {
private:
	// Extremalwerte aka Boundaries
	coordinate min, max;

	// Periodizität des Musters
	bool periodic;
public:
	nodelist();
	nodelist(bool periodicity);
	~nodelist();

	/* Listenspezifisches */
	// zeigt die Liste an.
	void display();

	/* Datenspezifisches */
	// Minimalwerte setzen
	void setMins(coordinate mins);
	// Maximalwerte setzen
	void setMaxs(coordinate maxs);
	// Bewegt die Liste um den gegebenen Vektor
	void shiftList(coordinate shifter);
	// Skaliert die Liste um die gegebenen Faktoren a,b, c
	void scaleList(double a);
	// Periodizität
	bool isPeriodic();
	// skaliert die Liste, sodass die Dichte der gewünschten entspricht
	void setDensity(double density);
	// berechnet die Punktdichte der Liste
	double getDensity();
	// berechnet das Volumen
	double getVolume();
	// überprüft ob ein Knoten an einem Punkt exisitiert
	node * getAt(coordinate point);
	// gibt die Minimalwerte der Box zurück
	coordinate getMins();
	// gibt die Maximalwerte der Box zurück
	coordinate getMaxs();
	// gibt die Größe der Box zurück
	coordinate getLengths();
	// gibt die größte Featuresize zurück, die Sinn macht (für periodische Randbedingungen)
	double getMaxFeatureSize();
	// Gibt die 26 Vektoren zurück um das Muster periodisch fortzusetzen. TODO: n-Dimensional
	vector<coordinate> getShifters();
	// Gibt die 26 Vektoren zurück um das Muster periodisch fortzusetzen. Es werden nur diejenigen Verschiebungen zurückgegeben, in denen die Box um den Mittelpunkt mid mit der Seitenlänge 2*halfExtend liegt. Geht davon aus, das das Muster den Schwerpunkt im Ursprung hat. TODO: n-Dimensional
	vector<coordinate> getShifted(coordinate mid, double halfExtend);
	// Zählt die Punkte in einer gegebenen Kugel
	int pointsInside(coordinate mid, double r);
	// Zählt die Punkte in einer gegebenen Kugel
	int pointsInsidePeriodic(coordinate mid, double r);
	// gibt Statistiken der Liste als String zurück TODO: was noch?
	string listStats(const char commentDelimeter[] = "\t");

	/* Berechnungen */
	void neighbourDistribution();
	void lengthDistribution();
	void angleDistribution();
	void hyperuniformity();
};

#endif /* NODELIST_HPP_ */
