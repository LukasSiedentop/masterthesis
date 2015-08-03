/*
 * nodelist.cpp
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

#include "nodelist.hpp"

NodeHead::NodeHead() :
		first(NULL), last(NULL), size(-1) {
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
		first = node;
		last = node;
	}


}

Node::Node() :
		head(NULL), prev(NULL), next(NULL), x(0), y(0), z(0) {
}

Node::Node(double x, double y, double z) :
		head (NULL), prev(NULL), next(NULL), x(x), y(y), z(z) {
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

void Node::setX(double newX){
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

double Node::getX(){
	return x;
}
