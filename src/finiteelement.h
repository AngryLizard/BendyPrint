
#pragma once

#include "element.h"

using Tensor = Eigen::Matrix3d;
using Vec3 = Eigen::RowVector3d;

class FiniteElement : public Element
{

public:

	FiniteElement(const Eigen::MatrixXd& V, const Eigen::VectorXi& vertices, const Eigen::VectorXi& faces, double shear, double bulk);
	virtual ~FiniteElement();
	
	/////////////////////////////////////////////////////
public:

	virtual void compute(const Eigen::MatrixXd& V) override;
	virtual double diffTest(const Eigen::MatrixXd& V, double h) override;

	double diffTestE(const Eigen::MatrixXd& x, double h);
	double diffTestF(const Eigen::MatrixXd& x, double h);
	double diffTestG(const Eigen::MatrixXd& x, double h);

	/////////////////////////////////////////////////////
protected:

	virtual Eigen::MatrixXd computeDeformation(const Eigen::MatrixXd& x) const = 0;
	virtual Eigen::MatrixXd computeStrain(const Eigen::MatrixXd& F) const = 0;

	virtual Eigen::MatrixXd dUdx(const Eigen::MatrixXd& dUdF) const = 0;
	virtual Eigen::MatrixXd dEdF(const Eigen::MatrixXd& F) const = 0;

	double _volume;

	/////////////////////////////////////////////////////
private:

	Eigen::MatrixXd dUdE(const Eigen::MatrixXd& E) const;
	double energy(const Eigen::MatrixXd& E) const;

	Eigen::MatrixXd _F, _E;
	double _shear, _bulk;
};