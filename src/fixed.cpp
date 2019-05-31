#include "fixed.h"

#include <Eigen/Eigen>


Fixed::Fixed(const Eigen::MatrixXd& V, const Eigen::VectorXi& indices, double k)
	: Element(V, indices, Eigen::VectorXi()), _k(k)
{
}

Fixed::~Fixed()
{

}


void Fixed::compute(const Eigen::MatrixXd& V, const Eigen::MatrixXd& ext)
{
	Eigen::MatrixXd x = extract(V);
	_energy = energy(x);
	_G = computeGradient(x);

#ifdef COMPUTE_FINITE_DIFFERENCE
	_fG = computeGradient(x);
#endif
}

double Fixed::diffTest(const Eigen::MatrixXd& V, double h)
{
	Eigen::MatrixXd x = extract(V);
	Eigen::MatrixXd G = computeGradient(x);
	const double e = energy(x);

	Eigen::MatrixXd hG = Eigen::MatrixXd::Zero(x.rows(), x.cols());
	for (int i = 0; i < x.rows(); i++)
	{
		for (int j = 0; j < x.cols(); j++)
		{
			Eigen::MatrixXd hX = x;
			hX(i, j) += h;

			hG(i, j) = (energy(hX) - e) / h;
		}
	}

	return (hG - G).squaredNorm();
}

/////////////////////////////////////////////////////

Eigen::MatrixXd Fixed::computeGradient(const Eigen::MatrixXd& x) const
{
	return _k * (x - _X);
}

double Fixed::energy(const Eigen::MatrixXd& x) const
{
	return 0.5 * _k * (x - _X).rowwise().squaredNorm().sum();
}