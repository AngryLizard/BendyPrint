
#pragma once

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <iostream>

#define EPS 0.0001f
//#define COMPUTE_FINITE_DIFFERENCE

class Element
{

public:

	Element(const Eigen::MatrixXd& V, const Eigen::VectorXi& vertices, const Eigen::VectorXi& faces);
	virtual ~Element();

	double getEnergy() const;

	/////////////////////////////////////////////////////
public:

	virtual void compute(const Eigen::MatrixXd& V) = 0;
	virtual double diffTest(const Eigen::MatrixXd& V, double h) = 0;

	void addGradient(Eigen::MatrixXd& G) const;
	void addHessian(Eigen::MatrixXd& H) const;

	void addFiniteGradient(Eigen::MatrixXd& G) const;

	void colorFaces(const Eigen::Vector3d& rgb, Eigen::MatrixXd& C) const;

	/////////////////////////////////////////////////////
protected:

	Eigen::MatrixXd extract(const Eigen::MatrixXd& V);

	Eigen::VectorXi _indices, _faces;
	Eigen::MatrixXd _X, _G, _fG, _H;
	double _energy;
};
