
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

	virtual Eigen::MatrixXd compute_F(const Eigen::MatrixXd& X) const override;
	virtual Eigen::MatrixXd compute_E(const Eigen::MatrixXd& F) const override;
	virtual double compute_U(const Eigen::MatrixXd& E) const override;
	
	virtual Eigen::MatrixXd compute_dUdE(const Eigen::MatrixXd& E) const override;
	virtual Eigen::MatrixXd compute_dEdF(const Eigen::MatrixXd& F) const override;
	virtual Eigen::MatrixXd compute_dFdx(const Eigen::MatrixXd& x, int32_t i, int32_t j) const override;

	virtual Eigen::MatrixXd compute_ddFddx(const Eigen::MatrixXd& x) const override;
	virtual Eigen::MatrixXd compute_ddEddF(const Eigen::MatrixXd& F) const override;
	virtual Tensor compute_ddUddE(const Eigen::MatrixXd& E) const override;


	/////////////////////////////////////////////////////
private:
	Rotation _invdX;
};
