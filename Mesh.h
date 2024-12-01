#pragma once
#include <iostream>
#include <vector>
#include <array>

class Element {
public:
	Element(double L, double E, double A, double k, double p_bar, double P = 0);

	double GetStiffness(size_t i, size_t j) const;
	double GetForce(size_t i) const;

protected:

	std::array<std::array<double, 2>, 2> stiffnessMatrix;
	std::array<double, 2> forceVector;
};

class Mesh
{
public:
	Mesh(double length, int numElements, double stiffness, double area, double spring, 
		double distributedLoad, double concentratedLoad);

	void RunAnalysis();

protected:

	void Discretize();
	void AssembleStiffnessMatrix();
	void AssembleForceVector();
	void ApplyBC();
	void Solve();
	void PrintResults() const;

	void PrintMatrix(const std::vector<std::vector<double>>& mat) const;
	void PrintVector(const std::vector<double>& vec) const;

	void SetStiffness(size_t row, size_t col, double val);
	void SetForce(size_t row, double val);

	double L;
	int n;
	double E;
	double A;
	double k;
	double p_bar;
	double P;

	std::vector<Element> elements;

	std::vector<std::array<int, 2>> connectivityMatrix;
	std::vector<std::array<double, 2>> coordinateMatrix;

	std::vector<std::vector<double>> globalStiffnessMatrix;
	std::vector<double> globalDisplacement;
	std::vector<double> globalForce;
};


