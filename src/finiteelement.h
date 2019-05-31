
#pragma once

#include "element.h"
#include "Tensor.h"

using Rotation = Eigen::Matrix3d;
using Vec3 = Eigen::RowVector3d;

class FiniteElement : public Element
{

public:

	FiniteElement(const Eigen::MatrixXd& V, const Eigen::VectorXi& vertices, const Eigen::VectorXi& faces, double shear, double bulk);
	virtual ~FiniteElement();
	
	/////////////////////////////////////////////////////
public:

	virtual void compute(const Eigen::MatrixXd& V, const Eigen::MatrixXd& ext) override;
	virtual double diffTest(const Eigen::MatrixXd& V, double h) override;

	// Gradient
	Eigen::MatrixXd diffTestdUdE(const Eigen::MatrixXd& E, double h) const;
	Eigen::MatrixXd diffTestdEdF(const Eigen::MatrixXd& F, double h) const;

	Tensor diffTestdFdx(const Eigen::MatrixXd& x, double h) const;
	Eigen::MatrixXd diffTestdUdx(const Eigen::MatrixXd& x, double h) const;
	Eigen::MatrixXd diffTestdUdF(const Eigen::MatrixXd& F, double h) const;

	// Hessian
	Eigen::MatrixXd diffTestddEddF(const Eigen::MatrixXd& F, double h) const;
	Tensor diffTestddUddE(const Eigen::MatrixXd& E, double h) const;

	Eigen::MatrixXd diffTestddUdEdF(const Eigen::MatrixXd& F, double h) const;
	Eigen::MatrixXd diffTestddUddF(const Eigen::MatrixXd& F, double h) const;
	Eigen::MatrixXd diffTestddUdFdx(const Eigen::MatrixXd& x, double h) const;
	Tensor diffTestddUdFdxT(const Eigen::MatrixXd& x, double h) const;
	Tensor diffTestddUddx(const Eigen::MatrixXd& x, double h) const;


	/////////////////////////////////////////////////////
protected:

	// Gradient
	virtual Eigen::MatrixXd compute_F(const Eigen::MatrixXd& x) const = 0;
	virtual Eigen::MatrixXd compute_E(const Eigen::MatrixXd& F) const = 0;
	virtual double compute_U(const Eigen::MatrixXd& E) const = 0;

	virtual Eigen::MatrixXd compute_dUdE(const Eigen::MatrixXd& E) const = 0;
	virtual Eigen::MatrixXd compute_dEdF(const Eigen::MatrixXd& F) const = 0;
	virtual Eigen::MatrixXd compute_dFdx(const Eigen::MatrixXd& x, int32_t i, int32_t j) const = 0;
	virtual Tensor compute_dFdx(const Eigen::MatrixXd& x) const;

	virtual Eigen::MatrixXd compute_dUdF(const Eigen::MatrixXd& dUdE, const Eigen::MatrixXd& dEdF) const;
	virtual Eigen::MatrixXd compute_dUdx(const Eigen::MatrixXd& dUdF, const Tensor& dFdx) const;

	// Hessian
	virtual Eigen::MatrixXd compute_ddFddx(const Eigen::MatrixXd& x) const = 0;
	virtual Eigen::MatrixXd compute_ddEddF(const Eigen::MatrixXd& F) const = 0;
	virtual Tensor compute_ddUddE(const Eigen::MatrixXd& E) const = 0;

	virtual Eigen::MatrixXd compute_ddUdEdF(const Eigen::MatrixXd& dEdF, const Tensor& ddUddE) const;
	virtual Eigen::MatrixXd compute_ddUddF(const Eigen::MatrixXd& dEdF, const Eigen::MatrixXd& ddUdEdF, const Eigen::MatrixXd& ddEddF, const Eigen::MatrixXd& dUdE) const;
	virtual Eigen::MatrixXd compute_ddUdFdx(const Tensor& dFdx, const Eigen::MatrixXd& ddUddF) const;
	virtual Tensor compute_ddUddx(const Tensor& dFdx, const Eigen::MatrixXd& ddUdFdx, const Eigen::MatrixXd& ddFddx, const Eigen::MatrixXd& dUdF) const;

	double _volume;
	double _shear, _bulk;

	/////////////////////////////////////////////////////
private:

	Eigen::MatrixXd _F, _E;
};
