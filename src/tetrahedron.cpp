#include "tetrahedron.h"

#include <Eigen/Eigen>


Tetrahedron::Tetrahedron(const Eigen::MatrixXd& V, const Eigen::VectorXi& indices, const Eigen::VectorXi& faces, double shear, double bulk)
: FiniteElement(V, indices, faces, shear, bulk)
{
	// Compute inverse initial deformation
	Rotation dX;
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


double Tetrahedron::compute_U(const Eigen::MatrixXd& E) const
{
	const double trE = E.trace();
	const double trSqr = E.cwiseProduct(E).sum();
	return (_shear / 2 * trE * trE + _bulk * trSqr) * _volume;
}

Eigen::MatrixXd Tetrahedron::compute_F(const Eigen::MatrixXd& x) const
{
	Rotation du;
	const auto& root = x.row(3);
	du << x.row(0) - root, x.row(1) - root, x.row(2) - root;
	return du * _invdX;
}

Eigen::MatrixXd Tetrahedron::compute_E(const Eigen::MatrixXd& F) const
{
	return (F.transpose() * F - Rotation::Identity()) / 2;
}


Eigen::MatrixXd Tetrahedron::compute_dUdE(const Eigen::MatrixXd& E) const
{
	Eigen::MatrixXd I = Rotation::Identity();
	return (_shear * E.trace() * I + 2 * _bulk * E) * _volume;
}

Eigen::MatrixXd Tetrahedron::compute_dEdF(const Eigen::MatrixXd& F) const
{
	return  F.transpose();
}

Eigen::MatrixXd Tetrahedron::compute_dFdx(const Eigen::MatrixXd& x, int32_t i, int32_t j) const
{
	Rotation dF = Rotation::Zero();
	if (i < 3)
	{
		dF(i, j) = 1.0;
	}
	else
	{
		dF.col(j) = Vec3(-1.0, -1.0, -1.0);
	}
	return dF * _invdX;
}

Eigen::MatrixXd Tetrahedron::compute_ddFddx(const Eigen::MatrixXd& x) const
{
	return _invdX;
}

Eigen::MatrixXd Tetrahedron::compute_ddEddF(const Eigen::MatrixXd& F) const
{
	return Eigen::Matrix3d::Identity();
}

Tensor Tetrahedron::compute_ddUddE(const Eigen::MatrixXd& E) const
{
	Tensor out(E.rows(), E.cols(), Eigen::Matrix3d::Zero());
	
	// Edge/Corner entries
	for (int i = 0; i < E.rows(); i++)
	{
		for (int j = 0; j < E.cols(); j++)
		{
			out(i, j)(i, j) += 2 * _bulk * _volume;
		}
	}

	// Diagonal entries
	for (int i = 0; i < E.rows(); i++)
	{
		for (int j = 0; j < E.cols(); j++)
		{
			out(i, i)(j, j) += _shear * _volume;
		}
	}

	return out;
}