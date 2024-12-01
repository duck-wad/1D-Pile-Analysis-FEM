#include <iomanip>
#include "Mesh.h"

Element::Element(double L, double E, double A, double k, double p_bar, double P) {

	double c1, c2 = 0;

	for (size_t i = 0; i <= 1; i++) {
		for (size_t j = 0; j <= 1; j++) {
			if (i == j) {
				c1 = 1.0;
				c2 = 1.0;
			}
			else {
				c1 = -1.0;
				c2 = 0.5;
			}
			stiffnessMatrix[i][j] = ((A * E) / L) * c1 + ((k * L) / 3) * c2;
		}

		if (i == 0) {
			forceVector[i] = -1 * (p_bar * L) / 2 + P;
		}
		else
			forceVector[i] = (p_bar * L) / 2 + P;
	}
}

double Element::GetStiffness(size_t i, size_t j) const {
	return stiffnessMatrix[i][j];
}

double Element::GetForce(size_t i) const {
	return forceVector[i];
}

Mesh::Mesh(double length, int numElements, double stiffness, double area, 
	double spring, double distributedLoad, double concentratedLoad)
	: L(length), n(numElements), E(stiffness), A(area), k(spring), 
	p_bar(distributedLoad), P(concentratedLoad) {

	Discretize();
}

void Mesh::Discretize() {
	double dx = L / n;
	elements.reserve(n);
	connectivityMatrix.reserve(n);
	coordinateMatrix.reserve(n);

	//for n elements there are n+1 nodes
	globalStiffnessMatrix.resize(n + 1, std::vector<double>(n + 1, 0.0));
	globalDisplacement.resize(n+1, 0.0);
	globalForce.resize(n+1, 0.0);

	double currentPos = 0.0;

	//populate connectivity matrix, coordinate matrix, element vector
	for (int i = 0; i < n; i++) {

		connectivityMatrix.push_back({ i, i + 1 });

		double startPos = currentPos;
		double endPos = startPos + dx;

		coordinateMatrix.push_back({ startPos, endPos });
		currentPos = endPos;

		//concentrated load is applied to top of pile which is the first node (zero index)
		if (i == 0) {
			elements.emplace_back(dx, E, A, k, p_bar, P);
		}
		else
			elements.emplace_back(dx, E, A, k, p_bar);
	}
}

void Mesh::RunAnalysis() {
	AssembleStiffnessMatrix();
	AssembleForceVector();
	ApplyBC();
	Solve();
	PrintResults();
}

void Mesh::AssembleStiffnessMatrix() {

	for (size_t i = 0; i < elements.size(); i++) {

		size_t node1 = connectivityMatrix[i][0];
		size_t node2 = connectivityMatrix[i][1];

		globalStiffnessMatrix[node1][node1] += elements[i].GetStiffness(0, 0);
		globalStiffnessMatrix[node1][node2] += elements[i].GetStiffness(0, 1);
		globalStiffnessMatrix[node2][node1] += elements[i].GetStiffness(1, 0);
		globalStiffnessMatrix[node2][node2] += elements[i].GetStiffness(1, 1);
	}
}

void Mesh::AssembleForceVector() {

	for (size_t i = 0; i < elements.size(); i++) {

		size_t node1 = connectivityMatrix[i][0];
		size_t node2 = connectivityMatrix[i][1];

		globalForce[node1] += elements[i].GetForce(0);
		globalForce[node2] += elements[i].GetForce(1);
	}
}

void Mesh::ApplyBC() {
	
	//assuming the bottom node is restrained with no movement
	double restraint = 0.0;
	size_t restrainedNode = static_cast<size_t>(n);

	//traverse row of restrainedNode of stiffness matrix, 
	//replacing all with zero except position of restrainedNode
	for (size_t i = 0; i < globalStiffnessMatrix.size(); i++) {
		if (i == restrainedNode) {
			SetStiffness(restrainedNode, i, 1.0);
		}
		else
			SetStiffness(restrainedNode, i, 0.0);
	}
	//set the restrained node force value to the prescribed displacement
	SetForce(restrainedNode, restraint);
}

void Mesh::SetStiffness(size_t row, size_t col, double val) {
	globalStiffnessMatrix[row][col] = val;
}

void Mesh::SetForce(size_t row, double val) {
	globalForce[row] = val;
}

void Mesh::Solve() {
	
	int size = globalForce.size();

	std::vector<double> lower(size - 1, 0.0);  // Lower diagonal
	std::vector<double> upper(size - 1, 0.0);  // Upper diagonal
	std::vector<double> diag(size, 0.0);       // Main diagonal

	// Extract diagonals from the global stiffness matrix
	for (size_t i = 0; i < size; i++) {
		diag[i] = globalStiffnessMatrix[i][i];  // Main diagonal
		if (i > 0) {
			lower[i - 1] = globalStiffnessMatrix[i][i - 1];  // Lower diagonal
		}
		if (i < size - 1) {
			upper[i] = globalStiffnessMatrix[i][i + 1];  // Upper diagonal
		}
	}

	// Forward sweep (elimination)
	for (size_t i = 1; i < size; i++) {
		double factor = lower[i - 1] / diag[i - 1];
		diag[i] -= factor * upper[i - 1];
		globalForce[i] -= factor * globalForce[i - 1];
	}

	// Back substitution
	globalDisplacement[size - 1] = globalForce[size - 1] / diag[size - 1];
	for (size_t i = size - 2; i != static_cast<size_t>(-1); i--) {
		globalDisplacement[i] = (globalForce[i] - upper[i] * globalDisplacement[i + 1]) / diag[i];
	}
}

void Mesh::PrintMatrix(const std::vector<std::vector<double>>& mat) const {

	for (size_t i = 0; i < mat.size(); ++i) {
		for (size_t j = 0; j < mat[i].size(); ++j) {
			std::cout << std::setw(10) << mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Mesh::PrintVector(const std::vector<double>& vec) const {
	for (size_t i = 0; i < vec.size(); ++i) {
		std::cout << "Node " << i << ": " << vec[i] << std::endl;
	}
	std::cout << std::endl;
}

void Mesh::PrintResults() const {

	/*std::cout << "Nodal displacements: " << std::endl;
	PrintVector(globalDisplacement);
	std::cout << std::endl;*/

	std::cout << "Surface settlement: " << globalDisplacement[0] << " [m]" << std::endl;
}