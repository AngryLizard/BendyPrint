#include "finiteelement.h"

#include <Eigen/Eigen>


FiniteElement::FiniteElement(const Eigen::MatrixXd& V, const Eigen::VectorXi& indices, const Eigen::VectorXi& faces, double shear, double bulk)
: Element(V, indices, faces), _volume(1.0), _shear(shear), _bulk(bulk)
{
}

FiniteElement::~FiniteElement()
{

}

void FiniteElement::compute(const Eigen::MatrixXd& V, const Eigen::MatrixXd& ext)
{
	const Eigen::MatrixXd x = extract(V);
	const Eigen::MatrixXd f = extract(ext);

	_F = compute_F(x);
	_E = compute_E(_F);
	_energy = compute_U(_E) - f.cwiseProduct(x - _X).sum();

	assert(_energy >= 0.0 && "Energy should always be positive.");

	const Tensor dFdx = compute_dFdx(x);
	const Eigen::MatrixXd dUdE = compute_dUdE(_E);
	const Eigen::MatrixXd dEdF = compute_dEdF(_F);
	const Eigen::MatrixXd dUdF = compute_dUdF(dUdE, dEdF);
	const Eigen::MatrixXd dUdx = compute_dUdx(dUdF, dFdx);
	_G = dUdx - f;

#ifdef COMPUTE_FINITE_DIFFERENCE
	_fG = diffTestG(x, 0.0001) - f;
#endif
}

double FiniteElement::diffTest(const Eigen::MatrixXd& V, double h)
{
	Eigen::MatrixXd x = extract(V);
	const Eigen::MatrixXd F = compute_F(x);
	const Eigen::MatrixXd E = compute_E(F);
	const double e = compute_U(E);

	// Compute first derivatives
	const Eigen::MatrixXd dUdE = compute_dUdE(E);
	const Eigen::MatrixXd dEdF = compute_dEdF(F);
	const Tensor dFdx = compute_dFdx(x);

	const Eigen::MatrixXd hdUdE = diffTestdUdE(E, h);
	const Eigen::MatrixXd hdEdF = diffTestdEdF(F, h);
	const Tensor hdFdx = diffTestdFdx(x, h);

	// Compute gradient
	const Eigen::MatrixXd dUdF = compute_dUdF(dUdE, dEdF);
	const Eigen::MatrixXd dUdx = compute_dUdx(dUdF, dFdx);

	const Eigen::MatrixXd hdUdF = diffTestdUdF(F, h);
	const Eigen::MatrixXd hdUdx = diffTestdUdx(x, h);

	// Compute second derivatives
	const Eigen::MatrixXd ddFddx = compute_ddFddx(x);
	const Eigen::MatrixXd ddEddF = compute_ddEddF(F);
	const Tensor ddUddE = compute_ddUddE(E);

	const Eigen::MatrixXd hddEddF = diffTestddEddF(F, h);
	const Tensor hddUddE = diffTestddUddE(E, h);

	// Compute Hessian
	const Eigen::MatrixXd ddUdEdF = compute_ddUdEdF(dEdF, ddUddE);
	const Eigen::MatrixXd ddUddF = compute_ddUddF(dEdF, ddUdEdF, ddEddF, dUdE);
	const Eigen::MatrixXd ddUdFdx = compute_ddUdFdx(dFdx, ddUddF);
	const Tensor ddUddx = compute_ddUddx(dFdx, ddUdFdx, ddFddx, dUdF);

	//const Tensor ddUdFdxT = diffTestddUdFdxT(x, h);

	const Eigen::MatrixXd hddUdEdF = diffTestddUdEdF(F, h);
	const Eigen::MatrixXd hddUddF = diffTestddUddF(F, h);
	const Eigen::MatrixXd hddUdFdx = diffTestddUdFdx(x, h);
	const Tensor hddUddx = diffTestddUddx(x, h);

	//std::cout << "...................." << std::endl;
	//std::cout << ddUddx.flatten() << std::endl;
	//std::cout << "- - - - - - - -" << std::endl;
	//std::cout << hddUddx.flatten() << std::endl;
	//std::cout << ddUddx.toString() << std::endl;
	//std::cout << "- - - - - - - -" << std::endl;
	//std::cout << hddUddx.toString() << std::endl;

	double err = 0.0;
	//err += (dFdx - hdFdx).squaredNorm();
	//err += (dEdF - hdEdF).squaredNorm();
	//err += (dUdE - hdUdE).squaredNorm();
	//err += (dUdF - hdUdF).squaredNorm();
	err += (dUdx - hdUdx).squaredNorm();
	//err += (ddEddF - hddEddF).squaredNorm();
	//err += (ddUddE - hddUddE).squaredNorm();
	//err += (ddUddF - hddUddF).squaredNorm();
	//err += (ddUdFdx - hddUdFdx).squaredNorm();
	//err += (ddUddx - hddUddx).squaredNorm();

	return err;
}

/////////////////////////////////////////////////////



Tensor FiniteElement::diffTestdFdx(const Eigen::MatrixXd& x, double h) const
{
	const Eigen::MatrixXd F = compute_F(x);

	Tensor fF(x.rows(), x.cols());
	for (int i = 0; i < x.rows(); i++)
	{
		for (int j = 0; j < x.cols(); j++)
		{
			Eigen::MatrixXd hx = x;
			hx(i, j) += h;

			const Rotation hF = compute_F(hx);
			fF[i][j] = (hF - F) / h;
		}
	}
	return fF;
}

Eigen::MatrixXd FiniteElement::diffTestdEdF(const Eigen::MatrixXd& F, double h) const
{
	const Eigen::MatrixXd E = compute_E(F);

	Tensor fE(F.rows(), F.cols());
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			Eigen::MatrixXd hF = F;
			hF(i, j) += h;

			const Rotation hE = compute_E(hF);
			fE[i][j] = (hE - E) / h;
		}
	}
	return fE * Eigen::MatrixXd::Identity(F.rows(), F.cols());
}

Eigen::MatrixXd FiniteElement::diffTestdUdE(const Eigen::MatrixXd& E, double h) const
{
	const double e = compute_U(E);

	Eigen::MatrixXd fU = Eigen::MatrixXd::Zero(E.rows(), E.cols());
	for (int i = 0; i < E.rows(); i++)
	{
		for (int j = 0; j < E.cols(); j++)
		{
			Eigen::MatrixXd hE = E;
			hE(i, j) += h;

			fU(i, j) = (compute_U(hE) - e) / h;
		}
	}
	return fU;
}

Eigen::MatrixXd FiniteElement::diffTestdUdF(const Eigen::MatrixXd& F, double h) const
{
	const Eigen::MatrixXd E = compute_E(F);
	const double e = compute_U(E);

	Eigen::MatrixXd fU = Eigen::MatrixXd::Zero(F.rows(), F.cols());
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			Eigen::MatrixXd hF = F;
			hF(i, j) += h;

			const Eigen::MatrixXd hE = compute_E(hF);
			fU(i, j) = (compute_U(hE) - e) / h;
		}
	}
	return fU;
}

Eigen::MatrixXd FiniteElement::diffTestdUdx(const Eigen::MatrixXd& x, double h) const
{
	const Eigen::MatrixXd F = compute_F(x);
	const Eigen::MatrixXd E = compute_E(F);
	const double e = compute_U(E);

	Eigen::MatrixXd fU = Eigen::MatrixXd::Zero(x.rows(), x.cols());
	for (int i = 0; i < x.rows(); i++)
	{
		for (int j = 0; j < x.cols(); j++)
		{
			Eigen::MatrixXd hx = x;
			hx(i, j) += h;

			const Eigen::MatrixXd hF = compute_F(hx);
			const Eigen::MatrixXd hE = compute_E(hF);
			fU(i, j) = (compute_U(hE) - e) / h;
		}
	}
	return fU;
}


Eigen::MatrixXd FiniteElement::diffTestddEddF(const Eigen::MatrixXd& F, double h) const
{
	const Eigen::MatrixXd dEdF = compute_dEdF(F);

	Tensor fE(F.rows(), F.cols());
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			Eigen::MatrixXd hF = F;
			hF(i, j) += h;

			const Rotation hdEdF = compute_dEdF(hF);
			fE[i][j] = (hdEdF - dEdF) / h;
		}
	}
	return fE * Eigen::MatrixXd::Identity(F.rows(), F.cols());
}

Tensor FiniteElement::diffTestddUddE(const Eigen::MatrixXd& E, double h) const
{
	const Eigen::MatrixXd dUdE = compute_dUdE(E);

	Tensor fU(E.rows(), E.cols());
	for (int i = 0; i < E.rows(); i++)
	{
		for (int j = 0; j < E.cols(); j++)
		{
			Eigen::MatrixXd hE = E;
			hE(i, j) += h;

			const Rotation hdUdE = compute_dUdE(hE);
			fU[i][j] = (hdUdE - dUdE) / h;
		}
	}
	return fU;
}

Eigen::MatrixXd FiniteElement::diffTestddUdEdF(const Eigen::MatrixXd& F, double h) const
{
	const Eigen::MatrixXd E = compute_E(F);
	const Eigen::MatrixXd dUdE = compute_dUdE(E);

	Tensor fU(F.rows(), F.cols());
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			Eigen::MatrixXd hF = F;
			hF(i, j) += h;

			const Eigen::MatrixXd hE = compute_E(hF);
			const Eigen::MatrixXd hdUdE = compute_dUdE(hE);
			fU(i, j) = (hdUdE - dUdE) / h;
		}
	}
	return fU * Eigen::MatrixXd::Identity(F.rows(), F.cols());
}

Eigen::MatrixXd FiniteElement::diffTestddUddF(const Eigen::MatrixXd& F, double h) const
{
	const Eigen::MatrixXd E = compute_E(F);
	const Eigen::MatrixXd dUdE = compute_dUdE(E);
	const Eigen::MatrixXd dEdF = compute_dEdF(F);
	const Eigen::MatrixXd dUdF = compute_dUdF(dUdE, dEdF);

	Tensor fU(F.rows(), F.cols());
	for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < F.cols(); j++)
		{
			Eigen::MatrixXd hF = F;
			hF(i, j) += h;

			const Eigen::MatrixXd hE = compute_E(hF);
			const Eigen::MatrixXd hdUdE = compute_dUdE(hE);
			const Eigen::MatrixXd hdEdF = compute_dEdF(hF);
			const Eigen::MatrixXd hdUdF = compute_dUdF(hdUdE, hdEdF);
			fU(i, j) = (hdUdF - dUdF) / h;
		}
	}
	return fU * Eigen::MatrixXd::Identity(F.rows(), F.cols());
}

Eigen::MatrixXd FiniteElement::diffTestddUdFdx(const Eigen::MatrixXd& x, double h) const
{
	const Eigen::MatrixXd F = compute_F(x);
	const Eigen::MatrixXd E = compute_E(F);
	const Eigen::MatrixXd dUdE = compute_dUdE(E);
	const Eigen::MatrixXd dEdF = compute_dEdF(F);
	const Eigen::MatrixXd dUdF = compute_dUdF(dUdE, dEdF);

	Tensor fU(x.rows(), x.cols());
	for (int i = 0; i < x.rows(); i++)
	{
		for (int j = 0; j < x.cols(); j++)
		{
			Eigen::MatrixXd hx = x;
			hx(i, j) += h;

			const Eigen::MatrixXd hF = compute_F(hx);
			const Eigen::MatrixXd hE = compute_E(hF);
			const Eigen::MatrixXd hdUdE = compute_dUdE(hE);
			const Eigen::MatrixXd hdEdF = compute_dEdF(hF);
			const Eigen::MatrixXd hdUdF = compute_dUdF(hdUdE, hdEdF);
			fU(i, j) = (hdUdF - dUdF) / h;
		}
	}
	return fU * Eigen::MatrixXd::Identity(F.rows(), F.cols());
}

Tensor FiniteElement::diffTestddUdFdxT(const Eigen::MatrixXd& x, double h) const
{
	const Eigen::MatrixXd F = compute_F(x);
	const Eigen::MatrixXd E = compute_E(F);
	const Eigen::MatrixXd dUdE = compute_dUdE(E);
	const Eigen::MatrixXd dEdF = compute_dEdF(F);
	const Eigen::MatrixXd dUdF = compute_dUdF(dUdE, dEdF);

	Tensor fU(x.rows(), x.cols());
	for (int i = 0; i < x.rows(); i++)
	{
		for (int j = 0; j < x.cols(); j++)
		{
			Eigen::MatrixXd hx = x;
			hx(i, j) += h;

			const Eigen::MatrixXd hF = compute_F(hx);
			const Eigen::MatrixXd hE = compute_E(hF);
			const Eigen::MatrixXd hdUdE = compute_dUdE(hE);
			const Eigen::MatrixXd hdEdF = compute_dEdF(hF);
			const Eigen::MatrixXd hdUdF = compute_dUdF(hdUdE, hdEdF);
			fU(i, j) = (hdUdF - dUdF) / h;
		}
	}
	return fU;
}

Tensor FiniteElement::diffTestddUddx(const Eigen::MatrixXd& x, double h) const
{
	const Eigen::MatrixXd F = compute_F(x);
	const Eigen::MatrixXd E = compute_E(F);
	const Eigen::MatrixXd dUdE = compute_dUdE(E);
	const Eigen::MatrixXd dEdF = compute_dEdF(F);
	const Tensor dFdx = compute_dFdx(x);
	const Eigen::MatrixXd dUdF = compute_dUdF(dUdE, dEdF);
	const Eigen::MatrixXd dUdx = compute_dUdx(dUdF, dFdx);

	Tensor fU(x.rows(), x.cols());
	for (int i = 0; i < x.rows(); i++)
	{
		for (int j = 0; j < x.cols(); j++)
		{
			Eigen::MatrixXd hx = x;
			hx(i, j) += h;

			const Eigen::MatrixXd hF = compute_F(hx);
			const Eigen::MatrixXd hE = compute_E(hF);
			const Eigen::MatrixXd hdUdE = compute_dUdE(hE);
			const Eigen::MatrixXd hdEdF = compute_dEdF(hF);
			const Tensor hdFdx = compute_dFdx(hx);
			const Eigen::MatrixXd hdUdF = compute_dUdF(hdUdE, hdEdF);
			const Eigen::MatrixXd hdUdx = compute_dUdx(hdUdF, hdFdx);
			fU(i, j) = (hdUdx - dUdx) / h;
		}
	}
	return fU;
}

/////////////////////////////////////////////////////


Tensor FiniteElement::compute_dFdx(const Eigen::MatrixXd& x) const
{
	Tensor out(x.rows(), x.cols());
	for (int i = 0; i < x.rows(); i++)
	{
		for (int j = 0; j < x.cols(); j++)
		{
			out[i][j] = compute_dFdx(x, i, j);
		}
	}
	return out;
}

Eigen::MatrixXd FiniteElement::compute_dUdx(const Eigen::MatrixXd& dUdF, const Tensor& dFdx) const
{
	return dFdx * dUdF;
}

Eigen::MatrixXd FiniteElement::compute_dUdF(const Eigen::MatrixXd& dUdE, const Eigen::MatrixXd& dEdF) const
{
	return dEdF.transpose() * dUdE;
}


Tensor FiniteElement::compute_ddUddx(const Tensor& dFdx, const Eigen::MatrixXd& ddUdFdx, const Eigen::MatrixXd& ddFddx, const Eigen::MatrixXd& dUdF) const
{
	return dFdx.mult(ddUdFdx);// +ddFddx.transpose() * dUdF; // dFdx was constant, second derivative should therefore be 0
}

Eigen::MatrixXd FiniteElement::compute_ddUdFdx(const Tensor& dFdx, const Eigen::MatrixXd& ddUddF) const
{
	return dFdx * ddUddF;
}

Eigen::MatrixXd FiniteElement::compute_ddUddF(const Eigen::MatrixXd& dEdF, const Eigen::MatrixXd& ddUdEdF, const Eigen::MatrixXd& ddEddF, const Eigen::MatrixXd& dUdE) const
{
	return dEdF.transpose() * ddUdEdF + ddEddF.transpose() * dUdE;
}

Eigen::MatrixXd FiniteElement::compute_ddUdEdF(const Eigen::MatrixXd& dEdF, const Tensor& ddUddE) const
{
	return ddUddE * dEdF.transpose();
}

/////////////////////////////////////////////////////