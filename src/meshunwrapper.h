
#pragma once

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <vector>

#include "tetrahedron.h"
#include "fixed.h"

#include "slicer.h"


class MeshUnwrapper
{
	/////////////////////////////////////////////////////
public:

	// Computes unwrapping from a given original and shell mesh given an altitude for the unwrapped mesh and a predicate for whether two faces can be connected (defaults to always)
	MeshUnwrapper(const Eigen::MatrixXd& oV, const Eigen::MatrixXi& oF, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double altitude, std::function<bool(int32_t, int32_t)> canConnect = [](int32_t, int32_t) {return true; });
	~MeshUnwrapper();

	/////////////////////////////////////////////////////
public:

	// Render unwrap
	void renderPlane(igl::opengl::glfw::Viewer& viewer, bool hasCeiling) const;
	
	// Module getters
	const Slicer* getSlicer() const;

	/////////////////////////////////////////////////////
private:

	Slicer* _slicer;

	Eigen::MatrixXd _pV;
	Eigen::MatrixXi _pF; // Plane

	Eigen::MatrixXd _dV;
	Eigen::MatrixXi _dF; // Ceiling
};
