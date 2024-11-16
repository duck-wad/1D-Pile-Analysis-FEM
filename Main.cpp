#include <iostream>
#include <vector>
#include "Mesh.h"

//This program is to simulate the settlement of a 1D pile under axial load using FEA
//Concentrated load and distributed load is applied to the pile
//Assume the load is applied to the top node and the bottom node is on a rigid rock layer with u=0

int main() {
	
	//Define pile parameters
	double pileLength = 20.0;	//m
	double stiffness = 1e5;		//Pa
	double area = 0.8;			//m^2
	double spring = 500;		//N/m

	//Define load [N]
	double distributedLoad = 100;			//N/m
	double concentratedLoad = 1000.0;		//N

	//Define number of elements
	int numElements = 100;

	//Initialize the mesh
	Mesh mesh(pileLength, numElements, stiffness, area, spring, distributedLoad, concentratedLoad);
	mesh.RunAnalysis();

	return 0;
}