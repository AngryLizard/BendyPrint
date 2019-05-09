#include "element.h"

#include <Eigen/Eigen>


Element::Element(const Eigen::MatrixXd& V, const Eigen::VectorXi& indices, const Eigen::VectorXi& faces)
: _indices(indices), _faces(faces)
{
	_X = extract(V);
}

Element::~Element()
{

}

double Element::getEnergy() const
{
	return _energy;
}

/////////////////////////////////////////////////////

void Element::addGradient(Eigen::MatrixXd& G) const
{
	const int32_t n = (int32_t)_indices.size();
	for (int32_t i = 0; i < n; i++)
	{
		G.row(_indices[i]) += _G.row(i);
	}
}

void Element::addHessian(Eigen::MatrixXd& H) const
{
	const int32_t n = (int32_t)_indices.size();
	for (int32_t i = 0; i < n; i++)
	{
		H.row(_indices[i]) += _H.row(i);
	}
}

void Element::colorFaces(const Eigen::Vector3d& rgb, Eigen::MatrixXd& C) const
{
	for (int i = 0; i < _faces.size(); i++)
	{
		C.row(_faces[i]) = rgb;
	}
}

/////////////////////////////////////////////////////

Eigen::MatrixXd Element::extract(const Eigen::MatrixXd& V)
{
	const int32_t n = (int32_t)_indices.size();
	Eigen::MatrixXd x = Eigen::MatrixXd::Zero(n, 3);
	for (int32_t i = 0; i < n; i++)
	{
		x.row(i) = V.row(_indices[i]);
	}
	return x;
}

