
#pragma once

#include <Eigen/Core>
#include <iostream>

#include "finiteelement.h"


class Tetrahedron : public FiniteElement
{

public:

	Tetrahedron(const Eigen::MatrixXd& V, const Eigen::VectorXi& indices, const Eigen::VectorXi& faces, double shear, double bulk);
	virtual ~Tetrahedron();

	/////////////////////////////////////////////////////
protected:

	virtual Eigen::MatrixXd computeDeformation(const Eigen::MatrixXd& X) const override;
	virtual Eigen::MatrixXd computeStrain(const Eigen::MatrixXd& F) const override;

	virtual Eigen::MatrixXd dUdx(const Eigen::MatrixXd& dUdF) const override;
	virtual Eigen::MatrixXd dUdF(const Eigen::MatrixXd& E, const Eigen::MatrixXd& F) const override;

	/////////////////////////////////////////////////////
private:
	Tensor _invdX;
};
