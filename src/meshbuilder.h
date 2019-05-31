
#pragma once

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <vector>

#include "tetrahedron.h"
#include "fixed.h"

#include "meshunwrapper.h"
#include "FEM.h"


class Meshbuilder
{
	/////////////////////////////////////////////////////
public:

	Meshbuilder(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	~Meshbuilder();

	// Generate shell from mesh
	void buildFromPlane(double depth, double shear, double bulk);

	// Assume shell already being generated, vertices need to be ordered correctly
	void buildFromVolume(double shear, double bulk);
	
	// Set FEM properties for faces
	void setExternalForce(int32_t face, const Eigen::Vector3d& force);
	void addFixedFace(int32_t face);
	void removeFixedFace(int32_t face);

	// Vertex location getter/setter
	void setVertexLocation(int32_t vertex, const Eigen::Vector3d& position);
	Eigen::Vector3d getVertexLocation(int32_t vertex) const;

	// Unwrap to plane
	void unwrap(double altitude);

	// Slice unwrapped
	void slice();

	// Run FEM simulation for a given amount of iterations
	void simulate(int32_t iterations);

	// Tests FEM simulation, results are printed to stdout
	void testSimulation();

	/////////////////////////////////////////////////////
public:

	// Render without coloration
	void renderShell(igl::opengl::glfw::Viewer& viewer) const;

	// Render with FEM coloration
	void renderFEM(igl::opengl::glfw::Viewer& viewer, bool constraints) const;

	// Module getters
	const FEM* getFEM() const;
	const MeshUnwrapper* getUnwrapper() const;

	/////////////////////////////////////////////////////
private:

	FEM* _fem;
	MeshUnwrapper* _unwrapper;

	double _depth;

	Eigen::MatrixXd _V;
	Eigen::MatrixXi _F;
	Eigen::MatrixXi _S; // Shell
	Eigen::MatrixXd _C; // Colour

	Eigen::MatrixXd _oV;
	Eigen::MatrixXi _oF; // Origin

	Eigen::MatrixXd _pV;
	Eigen::MatrixXi _pF; // Plane

	Eigen::MatrixXd _dV;
	Eigen::MatrixXi _dF; // Ceiling
};
