
#pragma once

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <vector>

#include "tetrahedron.h"
#include "fixed.h"

#include "IntervalSlice.h"

struct Face
{
	int32_t face; // Face to be projected
	int32_t faceEdge; // Face edge to attach
	int32_t prev; // Face to be attached to
	int32_t prevEdge; // Edge to be attach to
};

class Meshbuilder
{
	/////////////////////////////////////////////////////
public:

	Meshbuilder(const Eigen::MatrixXd& V,const Eigen::MatrixXi& F, double depth);
	Meshbuilder(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	~Meshbuilder();

	double totalEnergy() const;
	void addFixedVertex(int32_t vertex, const Eigen::Vector3d& position);
	void setFixedVertex(int32_t vertex, const Eigen::Vector3d& position);
	Eigen::Vector3d getFixedVertex(int32_t vertex) const;

	/////////////////////////////////////////////////////
public:

	void computeElements();
	void gradientTest(double h);
	void gradientDescent(double lambda, int32_t maxIterations);
	void renderShell(igl::opengl::glfw::Viewer& viewer) const;

	/////////////////////////////////////////////////////
public:

	void computePlane(bool hasCeiling, double altitude);
	void renderPlane(igl::opengl::glfw::Viewer& viewer) const;

	/////////////////////////////////////////////////////
public:

	IntervalBread createIntervalBread(const Eigen::Vector2d& dir, double density) const;
	bool lineIntersection(double y, const Eigen::Vector2d& a, const Eigen::Vector2d& b, Interval& interval) const;
	bool tracePlane(IntervalBread& bread, double density) const;

	void renderSlices(const IntervalBread& bread, igl::opengl::glfw::Viewer& viewer, double altitude) const;

	/////////////////////////////////////////////////////
private:

	std::vector<Element*> _elements;
	Eigen::MatrixXd _gradient;
#ifdef COMPUTE_FINITE_DIFFERENCE
	Eigen::MatrixXd _finiteGradient;
#endif

	double _energy;

	Eigen::MatrixXd _V;
	Eigen::MatrixXi _F;
	Eigen::MatrixXi _S; // Shell
	Eigen::MatrixXd _C; // Colour

	Eigen::MatrixXd _oV;
	Eigen::MatrixXi _oF; // Origin

	Eigen::MatrixXd _pV;
	Eigen::MatrixXi _pF; // Plane
};
