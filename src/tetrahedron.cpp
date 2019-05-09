#include "tetrahedron.h"

#include <Eigen/Eigen>


Tetrahedron::Tetrahedron(const Eigen::MatrixXd& V, const Eigen::VectorXi& indices, const Eigen::VectorXi& faces, double shear, double bulk)
: FiniteElement(V, indices, faces, shear, bulk)
{
	// Compute inverse initial deformation
	Tensor dX;
	const auto root = _X.row(3);
	dX << _X.row(0) - root, _X.row(1) - root, _X.row(2) - root;
	_invdX = dX.inverse();

	// Compute volume
	_volume = dX.determinant() / 6;
	if (_volume < 0.0)
	{
		std::cout << "Negative volume, invalid configuration" << std::endl;
	}
}

Tetrahedron::~Tetrahedron()
{

}


Eigen::MatrixXd Tetrahedron::computeDeformation(const Eigen::MatrixXd& x) const
{
	Tensor du;
	const auto& root = x.row(3);
	du << x.row(0) - root, x.row(1) - root, x.row(2) - root;
	return du * _invdX;
}

Eigen::MatrixXd Tetrahedron::computeStrain(const Eigen::MatrixXd& F) const
{
	return (F.transpose() * F - Tensor::Identity()) / 2;
}

Eigen::MatrixXd Tetrahedron::dUdx(const Eigen::MatrixXd& dUdF) const
{
	const Tensor Tinv = _invdX.transpose();
	const Tensor U = dUdF * Tinv;

	Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(4, 3);
	//gradient << f1, f2, f3, ;
	gradient.block<3,3>(0,0) = U;
	gradient.row(3) = -U.colwise().sum();

	return gradient;
}

Eigen::MatrixXd Tetrahedron::dEdF(const Eigen::MatrixXd& F) const
{
	return F;
}