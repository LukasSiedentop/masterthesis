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
	// Konstruktor für eine leere Liste mit gegebener Periodizität
	nodelist(bool periodicity);
	// Konstruiert eine nichtperiodische Nodelist aus einem Vektor von Nodes.
	nodelist(std::vector<node>& vec, bool periodicity);
	~nodelist();

	/* Listenspezifisches */
	// zeigt die Liste an.
	void display();

	/* Datenspezifisches */
	// Gitb die Liste als einfachen Vektor zurück mit den Daten der Knoten (keine Pointer)
	std::vector<class node> getVector();
	// Minimalwerte setzen
	void setMins(coordinate mins);
	// Maximalwerte setzen
	void setMaxs(coordinate maxs);
	// Bewegt die Liste um den gegebenen Vektor
	void shiftList(coordinate shifter);
	// Skaliert die Liste um die gegebenen Faktoren a,b, c
	void scaleList(double a);
	// Periodizität zurückgeben
	bool isPeriodic();
	// skaliert die Liste, sodass die Dichte der gewünschten entspricht
	void setDensity(double density);
	// berechnet die Punktdichte der Liste
	double getDensity();
	// berechnet das Volumen
	double getVolume();
	// gibt den Knoten an einem Punkt zurück sofern er exisitiert
	node * getAt(coordinate point);
	// gibt die Minimalwerte der Box zurück
	coordinate getMins();
	// gibt die Maximalwerte der Box zurück
	coordinate getMaxs();
	// gibt die Größe der Box zurück
	coordinate getLengths();
	// gibt die Mitte der Box zurück
	coordinate getMid();
	// gibt die größte Featuresize zurück, die Sinn macht (für periodische Randbedingungen)
	double getMaxFeatureSize();
	// Gibt die 26 Vektoren zurück um das Muster periodisch fortzusetzen. TODO: n-Dimensional
	vector<coordinate> getShifters();
	// Gibt die Vektoren zurück um das Muster periodisch fortzusetzen. Es werden nur diejenigen Verschiebungen zurückgegeben, in denen die Box um den Mittelpunkt mid mit der Seitenlänge 2*halfExtend liegt. TODO: n-Dimensional
	vector<coordinate> getShifted(coordinate mid, double halfExtend);
	// Zählt die Punkte in einer gegebenen Kugel
	int pointsInside(coordinate mid, double r, double rSqr);
	// Zählt die Punkte in einer gegebenen Kugel
	int pointsInsidePeriodic(coordinate mid, double r);
	// gibt Statistiken der Liste als String zurück
	string listStats(const char commentDelimeter[] = "\t");
	// macht das Muster vergleichbar mit anderen: Punktdichte=1, Mittelpunkt=(0,0,0)
	double normalize();

	/* Berechnungen */
	void neighbourDistribution();
	void lengthDistribution();
	void angleDistribution();
	void hyperuniformity();
	void writePOV();
};

#endif /* NODELIST_HPP_ */
