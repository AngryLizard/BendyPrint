#include "meshbuilder.h"
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/combine.h>

#include <set>
#include <queue>

Meshbuilder::Meshbuilder(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
	: _oV(V), _oF(F), _depth(0.0)
{
	_fem = new FEM();
	_unwrapper = nullptr;
}

Meshbuilder::~Meshbuilder()
{
	delete _fem;
	if (_unwrapper) delete _unwrapper;
}

void Meshbuilder::buildFromPlane(double depth, double shear, double bulk)
{
	_depth = depth;

	Eigen::MatrixXd vN;
	igl::per_vertex_normals(_oV, _oF, vN);

	// Make vertices
	const int32_t verts = _oV.rows();
	_V = Eigen::MatrixXd::Zero(verts * 2, 3);
	for (int i = 0; i < verts; i++)
	{
		_V.row(i) = _oV.row(i);
	}

	_fem->reset(verts * 2);

	// Make faces and elements
	const int32_t faces = _oF.rows();
	_F = Eigen::MatrixXi(faces * 9, 3);
	_C = Eigen::MatrixXd(faces * 9, 3);
	for (int i = 0; i < faces; i++)
	{
		// Compute face indices
		Eigen::Vector3i topf(i, i + faces * 2, i + faces * 3);
		Eigen::Vector3i celf(i + faces, i + faces * 4, i + faces * 5);
		Eigen::Vector3i midf(i + faces * 6, i + faces * 7, i + faces * 8);

		// compute floor and ceiling
		_F.row(topf.x()) = _oF.row(i);
		const Eigen::Vector3i& floor = _F.row(topf.x());
		_F.row(celf.x()) = floor + Eigen::Vector3i::Constant(verts);
		const Eigen::Vector3i& ceil = _F.row(celf.x());
		std::swap(_F.row(celf.x()).x(), _F.row(celf.x()).y());

		// Temporarily set to undeformed state by using face normal
		_V.row(ceil.x()) = _V.row(floor.x()) - vN.row(floor.x()) * depth;
		_V.row(ceil.y()) = _V.row(floor.y()) - vN.row(floor.y()) * depth;
		_V.row(ceil.z()) = _V.row(floor.z()) - vN.row(floor.z()) * depth;

		// Create elements
		_fem->addElement(new Tetrahedron(_V, Eigen::Vector4i(floor.x(), floor.y(), floor.z(), ceil.x()), topf, shear, bulk));
		_fem->addElement(new Tetrahedron(_V, Eigen::Vector4i(ceil.y(), ceil.x(), ceil.z(), floor.z()), celf, shear, bulk));
		_fem->addElement(new Tetrahedron(_V, Eigen::Vector4i(ceil.x(), ceil.y(), floor.y(), floor.z()), midf, shear, bulk));

		// Create visual walls
		_F.row(topf.y()) = Eigen::Vector3i(floor.x(), ceil.x(), floor.y());
		_F.row(topf.z()) = Eigen::Vector3i(floor.x(), floor.z(), ceil.x());

		_F.row(celf.y()) = Eigen::Vector3i(floor.y(), ceil.x(), ceil.y());
		_F.row(celf.z()) = Eigen::Vector3i(floor.y(), ceil.y(), ceil.z());

		_F.row(midf.x()) = Eigen::Vector3i(floor.z(), floor.y(), ceil.z());
		_F.row(midf.y()) = Eigen::Vector3i(floor.z(), ceil.z(), ceil.x());
		_F.row(midf.z()) = Eigen::Vector3i(floor.y(), floor.z(), ceil.x());
	}
}

void Meshbuilder::buildFromVolume(double shear, double bulk)
{
	_depth = 0.0;

	// Make vertices
	const int32_t verts = _oV.rows();
	_V = Eigen::MatrixXd::Zero(verts, 3);
	for (int i = 0; i < verts; i++)
	{
		_V.row(i) = _oV.row(i);
	}
	_fem->reset(verts);

	// Make elements
	_F = Eigen::MatrixXi(verts, 3);
	_C = Eigen::MatrixXd(verts, 3);
	for (int i = 0; i < verts; i += 4)
	{
		const Eigen::Vector4i idx = Eigen::Vector4i(i, i + 1, i + 2, i + 3);

		_F.row(i) = Eigen::Vector3i(i, i + 1, i + 2);
		_F.row(i + 1) = Eigen::Vector3i(i + 1, i + 3, i + 2);
		_F.row(i + 2) = Eigen::Vector3i(i, i + 1, i + 3);
		_F.row(i + 3) = Eigen::Vector3i(i, i + 3, i + 2);

		// Create elements
		_fem->addElement(new Tetrahedron(_V, idx, idx, shear, bulk));
	}
}


void Meshbuilder::setExternalForce(int32_t face, const Eigen::Vector3d& force)
{
	const Eigen::Vector3i& indices = _F.row(face);
	_fem->setExternalForce(indices.x(), force);
	_fem->setExternalForce(indices.y(), force);
	_fem->setExternalForce(indices.z(), force);
}

void Meshbuilder::addFixedFace(int32_t face)
{
	const Eigen::Vector3i& f = _F.row(face);

	_fem->addFixed(_V, f.x());
	_fem->addFixed(_V, f.y());
	_fem->addFixed(_V, f.z());
}

void Meshbuilder::removeFixedFace(int32_t face)
{
	const Eigen::Vector3i& f = _F.row(face);

	_fem->removeFixed(_V, f.x());
	_fem->removeFixed(_V, f.y());
	_fem->removeFixed(_V, f.z());
}


void Meshbuilder::setVertexLocation(int32_t vertex, const Eigen::Vector3d& position)
{
	_V.row(vertex) = position;
}

Eigen::Vector3d Meshbuilder::getVertexLocation(int32_t vertex) const
{
	return _V.row(vertex);
}

void Meshbuilder::unwrap(double altitude)
{
	// Update original vertex locations in case something changed
	for (int i = 0; i < _oF.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			const int32_t face = _F(i, j);
			_oV.row(face) = _V.row(face);
		}
	}

	if (_unwrapper) delete _unwrapper;
	_unwrapper = new MeshUnwrapper(_oV, _oF, _V, _F, altitude);
}

void Meshbuilder::simulate(int32_t iterations)
{
	_fem->gradientDescent(_V, 0.0001, iterations);
}

void Meshbuilder::testSimulation()
{
	_fem->gradientTest(_V, 0.00001);
}

/////////////////////////////////////////////////////

void Meshbuilder::renderShell(igl::opengl::glfw::Viewer& viewer) const
{
	viewer.data().clear();
	viewer.data().set_mesh(_V, _F);
}

void Meshbuilder::renderFEM(igl::opengl::glfw::Viewer& viewer, bool constraints) const
{
	_fem->computeElements(_V);
	if (constraints)
	{
		_fem->renderConstraints(viewer, _V, _F);
	}
	else
	{
		_fem->renderEnergies(viewer, _V, _F);
	}
}

const FEM* Meshbuilder::getFEM() const
{
	return _fem;
}

const MeshUnwrapper* Meshbuilder::getUnwrapper() const
{
	return _unwrapper;
}


/////////////////////////////////////////////////////