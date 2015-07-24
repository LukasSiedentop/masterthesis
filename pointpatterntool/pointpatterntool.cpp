/*
 * pointpatterntool.c
 *
 *  Created on: 24.07.2015
 *      Author: lukas
 */

#include <iostream>

int main(int argc, char *argv[]) {
	std::cout << "Hello World!" << std::endl;
	std::cout << "Given Arguments:" << std::endl;

	for (int i=0; i < argc; i++){
		std::cout << i << ": " << argv[i] << std::endl;
	}

	return 0;
}
