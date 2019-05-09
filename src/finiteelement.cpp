#include "finiteelement.h"

#include <Eigen/Eigen>


FiniteElement::FiniteElement(const Eigen::MatrixXd& V, const Eigen::VectorXi& indices, const Eigen::VectorXi& faces, double shear, double bulk)
: Element(V, indices, faces), _volume(1.0), _shear(shear), _bulk(bulk)
{
}

FiniteElement::~FiniteElement()
{

}

void FiniteElement::compute(const Eigen::MatrixXd& V)
{
	const Eigen::MatrixXd x = extract(V);
	_F = computeDeformation(x);
	_E = computeStrain(_F);
	_energy = energy(_E);

	const Eigen::MatrixXd dUdF = dUdE(_E) * dEdF(_F);
	_G = dUdx(dUdF);
}

double FiniteElement::diffTest(const Eigen::MatrixXd& V, double h)
{
	Eigen::MatrixXd x = extract(V);
	double err = 0.0;
	//err += diffTestE(x, h);
	err += diffTestF(x, h);
	//err += diffTestG(x, h);
	return err;
}

/////////////////////////////////////////////////////

double FiniteElement::diffTestE(const Eigen::MatrixXd& x, double h)
{
	const Eigen::MatrixXd F = computeDeformation(x);
	const Eigen::MatrixXd E = computeStrain(F);
	const Eigen::MatrixXd dU = dUdE(E);
	const double e = energy(E);

	Eigen::MatrixXd hU = Eigen::MatrixXd::Zero(E.rows(), E.cols());
	for (int i = 0; i < E.rows(); i++)
	{
		for (int j = 0; j < E.cols(); j++)
		{
			Eigen::MatrixXd hE = E;
			hE(i, j) += h;
			hU(i, j) = (energy(hE) - e) / h;
		}
	}
	return (hU - dU).squaredNorm();
}

double FiniteElement::diffTestF(const Eigen::MatrixXd& x, double h)
{
	const Eigen::MatrixXd F = computeDeformation(x);
	const Eigen::MatrixXd E = computeStrain(F);
	const Eigen::MatrixXd dU = dUdE(E) * dEdF(F);
	const double e = energy(E);

	Eigen::MatrixXd hU = Eigen::MatrixXd::Zero(F.rows(), F.cols());
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			Eigen::MatrixXd hF = F;
			hF(i, j) += h;

			const Eigen::MatrixXd hE = computeStrain(hF);
			hU(i, j) = (energy(hE) - e) / h;
		}
	}
	
	return (hU - dU).squaredNorm();
}

double FiniteElement::diffTestG(const Eigen::MatrixXd& x, double h)
{
	const Eigen::MatrixXd F = computeDeformation(x);
	const Eigen::MatrixXd E = computeStrain(F);
	const Eigen::MatrixXd dUdF = dUdE(E) * dEdF(F);
	const Eigen::MatrixXd dU = dUdx(dUdF);
	const double e = energy(E);

	Eigen::MatrixXd hU = Eigen::MatrixXd::Zero(x.rows(), x.cols());
	for (int i = 0; i < x.rows(); i++)
	{
		for (int j = 0; j < x.cols(); j++)
		{
			Eigen::MatrixXd hx = x;
			hx(i, j) += h;

			const Eigen::MatrixXd hF = computeDeformation(hx);
			const Eigen::MatrixXd hE = computeStrain(hF);
			hU(i, j) = (energy(hE) - e) / h;
		}
	}

	return (hU - dU).squaredNorm();
}

/////////////////////////////////////////////////////


Eigen::MatrixXd FiniteElement::dUdE(const Eigen::MatrixXd& E) const
{
	const double trE = E.trace();
	return (_shear * E.trace() * Tensor::Identity() + (2 * _bulk) * E.transpose()) * _volume;
}

double FiniteElement::energy(const Eigen::MatrixXd& E) const
{
	const double trE = E.trace();
	const double trSqr = E.cwiseProduct(E).sum();
	return (_shear / 2 * trE * trE + _bulk * trSqr) * _volume;
}
