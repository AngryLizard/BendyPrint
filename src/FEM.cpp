#include "FEM.h"
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/combine.h>

#include <set>
#include <queue>


FEM::FEM()
{
	_fixed = new Fixed(Eigen::MatrixXd(), Eigen::VectorXi(), 10.0);
	_elements.push_back(_fixed);
}

FEM::~FEM()
{
	// This also frees _fixed
	for (Element* element : _elements)
	{
		delete(element);
	}
}


void FEM::reset(int32_t n)
{
	_elements.clear();
	_elements.push_back(_fixed);

	if (_ext.rows() != n)
	{
		_ext = Eigen::MatrixXd::Zero(n, 3);
	}
}

void FEM::addElement(Element* element)
{
	_elements.push_back(element);
}

double FEM::totalEnergy() const
{
	return _energy;
}

double FEM::getEnergy(int32_t index) const
{
	return _elements[index]->getEnergy();
}


void FEM::setExternalForce(int32_t index, const Eigen::Vector3d& force)
{
	_ext.row(index) = force;
}

void FEM::addFixed(const Eigen::MatrixXd& V, int32_t index)
{
	std::vector<int32_t> indices;
	const Eigen::VectorXi Vi = _fixed->getVertices();
	for (int i = 0; i < Vi.size(); i++)
	{
		indices.push_back(Vi[i]);
	}

	if (std::find(indices.begin(), indices.end(), index) == indices.end()) indices.push_back(index);

	Eigen::VectorXi Wi(indices.size());
	for (int i = 0; i < indices.size(); i++)
	{
		Wi[i] = indices[i];
	}

	(*_fixed) = Fixed(V, Wi, 10000.0);
}

void FEM::removeFixed(const Eigen::MatrixXd& V, int32_t index)
{
	std::vector<int32_t> indices;
	const Eigen::VectorXi Vi = _fixed->getVertices();
	for (int i = 0; i < Vi.size(); i++)
	{
		indices.push_back(Vi[i]);
	}

	std::vector<int32_t>::iterator it;
	if ((it = std::find(indices.begin(), indices.end(), index)) != indices.end()) indices.erase(it);

	Eigen::VectorXi Wi(indices.size());
	for (int i = 0; i < indices.size(); i++)
	{
		Wi[i] = indices[i];
	}

	(*_fixed) = Fixed(V, Wi, 10000.0);
}

void FEM::computeElements(const Eigen::MatrixXd& V)
{
	// Compute element states
	const int32_t n = _elements.size();
	Eigen::VectorXd energies = Eigen::VectorXd::Zero(n);

	_gradient = Eigen::MatrixXd::Zero(V.rows(), 3);
	for (int32_t i = 0; i < n; i++)
	{
		Element* element = _elements[i];

		element->compute(V, _ext);
		element->addGradient(_gradient);
		energies[i] = element->getEnergy();

		assert(energies[i] >= 0.0 && "Energy should always be positive.");
	}

#ifdef COMPUTE_FINITE_DIFFERENCE
	_finiteGradient = Eigen::MatrixXd::Zero(_V.rows(), 3);
	for (int32_t i = 0; i < n; i++)
	{
		Element* element = _elements[i];
		element->addFiniteGradient(_finiteGradient);
	}

	/*
	std::cout << "__________________________" << std::endl;
	std::cout << _gradient << std::endl;
	std::cout << _finiteGradient << std::endl;
	*/
#endif

	_energy = energies.sum();
}

void FEM::gradientTest(const Eigen::MatrixXd& V, double h) const
{
	double err = 0.0;
	for (Element* element : _elements)
	{
		err += element->diffTest(V, h);
	}

	std::cout << "Tested! Error: " << err << std::endl;
}

void FEM::gradientDescent(Eigen::MatrixXd& V, double lambda, int32_t maxIterations)
{
	computeElements(V);
	double step = lambda;
	
	for (int i = 0; i < maxIterations; i++)
	{
		const auto oV = V;
		const auto oG = _gradient;

		//std::cout << ".........." << std::endl;
		//std::cout << oG << std::endl;

		const double energy = _energy;
		for (;;)
		{
			V = oV - oG * step;
			computeElements(V);

			if (_energy < energy)
			{
				step *= 1.1;
				break;
			}

			step *= 0.5;

			if (step < 1.0e-12)
			{
				V = oV;
				return;
			}
		}
	}
}

Eigen::MatrixXd FEM::computeColour(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) const
{
	// Get current energies
	const int32_t n = _elements.size();
	Eigen::MatrixXd C = Eigen::MatrixXd(F.rows(), 3);
	Eigen::VectorXd energies = Eigen::VectorXd::Zero(n);
	for (int32_t i = 0; i < n; i++)
	{
		energies[i] = _elements[i]->getEnergy();
	}

	// Compute colours
	const Eigen::Vector3d lo(1.0, 1.0, 1.0);
	const Eigen::Vector3d hi(1.0, 0.0, 0.0);

	const double min = energies.minCoeff();
	const double max = energies.maxCoeff();
	const double delta = max - min;
	if (delta > EPS)
	{
		for (int32_t i = 0; i < _elements.size(); i++)
		{
			const Element* element = _elements[i];
			const double alpha = (energies[i] - min) / delta;
			element->colorFaces(lo + alpha * (hi - lo), C);
		}
	}
	else
	{
		for (int32_t i = 0; i < _elements.size(); i++)
		{
			const Element* element = _elements[i];
			element->colorFaces(lo, C);
		}
	}
	return C;
}

void FEM::renderEnergies(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) const
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	viewer.data().set_face_based(false);

	Eigen::MatrixXd C = computeColour(V, F);
	viewer.data().set_colors(C);

	viewer.data().add_edges(V, V + _gradient, Eigen::RowVector3d(1.0, 1.0, 1.0));
}

void FEM::renderConstraints(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) const
{
	viewer.data().clear();
	viewer.data().set_mesh(V, F);

	viewer.data().set_face_based(false);

	// Colour external forces
	Eigen::MatrixXd C = Eigen::MatrixXd(V.rows(), 3);
	for (int i = 0; i < V.rows(); i++)
	{
		const Eigen::RowVectorXd& force = _ext.row(i);
		const double norm = force.norm();
		if (norm > EPS)
		{
			C.row(i) = force.cwiseAbs() / norm;
		}
		else
		{
			C.row(i) = Eigen::RowVector3d(1.0, 1.0, 1.0);
		}
	}

	// Colour fixed vertices blue
	const Eigen::VectorXi& vertices = _fixed->getVertices();
	for (int32_t i = 0; i < vertices.size(); i++)
	{
		C.row(vertices[i]) = Eigen::RowVector3d(0.0, 0.0, 0.0);
	}
	viewer.data().set_colors(C);
}
