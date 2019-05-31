
#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>

#include "bread.h"

#define EPS 0.0001

class Slicer
{
public:
	Slicer(const Eigen::MatrixXd& planeV, const Eigen::MatrixXi& planeF, const Eigen::MatrixXd& fillV, const Eigen::MatrixXi& fillF);
	virtual ~Slicer();

	// Convert 3D coordinates to 2D
	Eigen::Vector2d flatten(const Eigen::Vector3d& vec) const;

	// Convert 2D coordinates to 3D
	Eigen::Vector3d project(const Eigen::Vector2d& vec, double altitude) const;

	// Intervals
	/////////////////////////////////////////////////////

	// Create intervals from faces
	IntervalBread createIntervalBread(const Eigen::Vector2d& dir, double margin) const;

	// Intersect a line from a to b, min/max on given interval
	bool lineIntersection(const Eigen::Vector2d& a, const Eigen::Vector2d& b, Interval& interval) const;

	// Trace all faces on current position of IntervalBread
	bool tracePlane(IntervalBread& bread, double margin) const;

	// Create intervals from outlines
	IntervalBread createIntervalBread(const Eigen::Vector2d& dir, const OutlineBread& outlines, double margin) const;

	// Trace all outlines on current position of IntervalBread
	bool tracePlane(IntervalBread& bread, const OutlineBread& outlines, double margin) const;

	// Outlines
	/////////////////////////////////////////////////////

	// Creates outlines for a given plane
	OutlineBread createOutlineBread() const;

	// Iterates through all outline pieces while also implementing a margin. This margin makes it so additional lines are drawn around corners.
	void outlineSearch(const OutlineBread& bread, std::function<void(const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)> func, double margin) const;

	// Triangles
	/////////////////////////////////////////////////////

	// Create triangles for the fill
	TriangleBread createTriangleBread(double scale, double height) const;

	// Iterates through all triangle lines with given wall layers
	void triangleFillSearch(const TriangleBread& bread, std::function<void(const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)> func, int32_t layers, double margin) const;

	/////////////////////////////////////////////////////
private:

	Eigen::MatrixXd _planeV;
	Eigen::MatrixXi _planeF;
	Eigen::MatrixXd _fillV;
	Eigen::MatrixXi _fillF;
};