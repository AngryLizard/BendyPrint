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

	_G = dUdx(dUdF(_E, _F));

#ifdef COMPUTE_FINITE_DIFFERENCE
	_fG = diffTestG(x, 0.0001);
#endif
}

double FiniteElement::diffTest(const Eigen::MatrixXd& V, double h)
{
	Eigen::MatrixXd x = extract(V);
	const Eigen::MatrixXd F = computeDeformation(x);
	const Eigen::MatrixXd E = computeStrain(F);
	const double e = energy(E);

	const Eigen::MatrixXd ddUdF = dUdF(E, F);
	const Eigen::MatrixXd ddUdx = dUdx(ddUdF);

	//const Eigen::MatrixXd hdUdF = diffTestF(x, h);
	const Eigen::MatrixXd hdUdx = diffTestG(x, h);

	double err = 0.0;
	//err += (ddUdF - hdUdF).squaredNorm();
	err += (ddUdx - hdUdx).squaredNorm();

	return err;
}

/////////////////////////////////////////////////////

Eigen::MatrixXd FiniteElement::diffTestF(const Eigen::MatrixXd& x, double h)
{
	const Eigen::MatrixXd F = computeDeformation(x);
	const Eigen::MatrixXd E = computeStrain(F);
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
	
	return hU;
}

Eigen::MatrixXd FiniteElement::diffTestG(const Eigen::MatrixXd& x, double h)
{
	const Eigen::MatrixXd F = computeDeformation(x);
	const Eigen::MatrixXd E = computeStrain(F);
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

	return hU;
}

/////////////////////////////////////////////////////

double FiniteElement::energy(const Eigen::MatrixXd& E) const
{
	const double trE = E.trace();
	const double trSqr = E.cwiseProduct(E).sum();
	return (_shear / 2 * trE * trE + _bulk * trSqr) * _volume;
}
