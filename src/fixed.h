
#pragma once

#include <Eigen/Core>
#include <iostream>

#include "element.h"


class Fixed : public Element
{

public:

	Fixed(const Eigen::MatrixXd& V, const Eigen::VectorXi& indices, double k);
	virtual ~Fixed();
	
	virtual void compute(const Eigen::MatrixXd& V, const Eigen::MatrixXd& ext) override;
	virtual double diffTest(const Eigen::MatrixXd& V, double h) override;
	
	/////////////////////////////////////////////////////
private:

	Eigen::MatrixXd computeGradient(const Eigen::MatrixXd& x) const;
	double energy(const Eigen::MatrixXd& x) const;

	double _k;
};
