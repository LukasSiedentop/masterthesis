/*
 * pointpatterntool.c
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

#include <iostream>
#include "nodelist.hpp"

int main(int argc, char *argv[]) {
	std::cout << "Hello World!" << std::endl;
	std::cout << "Given Arguments:" << std::endl;

	for (int i=0; i < argc; i++){
		std::cout << i << ": " << argv[i] << std::endl;
	}


	Node * node = new Node(1, 2, 3);
	std::cout << "Node X: " << node->getX() << std::endl;

	NodeHead * list = new NodeHead();
	std::cout << "Listsize: " << list->length() << std::endl;
	list->add(node);
	std::cout << "Listsize: " << list->length() << std::endl;



	return 0;
}
