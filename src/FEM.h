
#pragma once

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <vector>

#include "tetrahedron.h"
#include "fixed.h"

#include "Slicer.h"


class FEM
{
	/////////////////////////////////////////////////////
public:

	FEM();
	~FEM();

	// Remove all elements but fixed
	void reset(int32_t n);

	// Add element
	void addElement(Element* element);

	// Compute energy
	double totalEnergy() const;
	double getEnergy(int32_t index) const;

	// Set external force for a vertex
	void setExternalForce(int32_t index, const Eigen::Vector3d& force);

	// Add/Remove fixed vertices
	void addFixed(const Eigen::MatrixXd& V, int32_t index);
	void removeFixed(const Eigen::MatrixXd& V, int32_t index);


	/////////////////////////////////////////////////////
public:

	// Compute gradient and energy for given vertex positions
	void computeElements(const Eigen::MatrixXd& V);

	// Test gradient computation
	void gradientTest(const Eigen::MatrixXd& V, double h) const;

	// Call gradient descent on vertex positions
	void gradientDescent(Eigen::MatrixXd& V, double lambda, int32_t maxIterations);

	// Compute colour depending on energy
	Eigen::MatrixXd computeColour(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) const;

	// Render relative energies
	void renderEnergies(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) const;

	// Render fixed elements and external forces
	void renderConstraints(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) const;

	/////////////////////////////////////////////////////
public:

	std::vector<Element*> _elements;
	Fixed* _fixed;

	Eigen::MatrixXd _ext;
	Eigen::MatrixXd _gradient;

#ifdef COMPUTE_FINITE_DIFFERENCE
	Eigen::MatrixXd _finiteGradient;
#endif

	double _energy;

};
