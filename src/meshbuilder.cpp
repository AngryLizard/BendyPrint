#include "meshbuilder.h"
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/combine.h>

#include <set>
#include <queue>

Meshbuilder::Meshbuilder(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double depth)
{
	_oV = V;
	_oF = F;

	const double shear = 50.0;
	const double bulk = 50.0;

	Eigen::MatrixXd vN, fN;
	igl::per_vertex_normals(_oV, _oF, vN);
	igl::per_face_normals(_oV, _oF, fN);

	// Make vertices
	const int32_t verts = _oV.rows();
	_V = Eigen::MatrixXd::Zero(verts * 2, 3);
	for (int i = 0; i < verts; i++)
	{
		_V.row(i) = _oV.row(i);
	}

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
		const auto delta = fN.row(i) * depth;
		_V.row(ceil.x()) = _V.row(floor.x()) - delta;
		_V.row(ceil.y()) = _V.row(floor.y()) - delta;
		_V.row(ceil.z()) = _V.row(floor.z()) - delta;

		// Create elements
		_elements.push_back(new Tetrahedron(_V, Eigen::Vector4i(floor.x(), floor.y(), floor.z(), ceil.x()), topf, shear, bulk));
		_elements.push_back(new Tetrahedron(_V, Eigen::Vector4i(ceil.y(), ceil.x(), ceil.z(), floor.z()), celf, shear, bulk));
		_elements.push_back(new Tetrahedron(_V, Eigen::Vector4i(ceil.x(), ceil.y(), floor.y(), floor.z()), midf, shear, bulk));

		// Create visual walls
		_F.row(topf.y()) = Eigen::Vector3i(floor.x(), ceil.x(), floor.y());
		_F.row(topf.z()) = Eigen::Vector3i(floor.x(), floor.z(), ceil.x());
		
		_F.row(celf.y()) = Eigen::Vector3i(floor.y(), ceil.x(), ceil.y());
		_F.row(celf.z()) = Eigen::Vector3i(floor.y(), ceil.y(), ceil.z());

		_F.row(midf.x()) = Eigen::Vector3i(floor.z(), floor.y(), ceil.z());
		_F.row(midf.y()) = Eigen::Vector3i(floor.z(), ceil.z(), ceil.x());
		_F.row(midf.z()) = Eigen::Vector3i(floor.y(), floor.z(), ceil.x());
	}

	// Properly locate vertices
	for (int i = 0; i < verts; i++)
	{
		_V.row(i + verts) = _oV.row(i) - vN.row(i) * depth;
	}
}

Meshbuilder::Meshbuilder(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
	_oV = V;
	_oF = F;

	const double shear = 50.0;
	const double bulk = 50.0;
	
	// Make vertices
	const int32_t verts = _oV.rows();
	_V = Eigen::MatrixXd::Zero(verts, 3);
	for (int i = 0; i < verts; i++)
	{
		_V.row(i) = _oV.row(i);
	}

	// Make elements
	_F = Eigen::MatrixXi(verts, 3);
	_C = Eigen::MatrixXd(verts, 3);
	for (int i = 0; i < verts; i+=4)
	{
		const Eigen::Vector4i idx = Eigen::Vector4i(i, i + 1, i + 2, i + 3);

		_F.row(i) =		Eigen::Vector3i(i, i+1, i+2);
		_F.row(i+1) =	Eigen::Vector3i(i+1, i+3, i+2);
		_F.row(i+2) =	Eigen::Vector3i(i, i+1, i+3);
		_F.row(i+3) =	Eigen::Vector3i(i, i+3, i+2);

		// Create elements
		_elements.push_back(new Tetrahedron(_V, idx, idx, shear, bulk));
	}

}

Meshbuilder::~Meshbuilder()
{
	for (Element* element : _elements)
	{
		delete(element);
	}
}

double Meshbuilder::totalEnergy() const
{
	return _energy;
}

void Meshbuilder::addFixedVertex(int32_t vertex, const Eigen::Vector3d& position)
{
	setFixedVertex(vertex, position);
	_elements.push_back(new Fixed(_V, Eigen::VectorXi::Constant(1, vertex), 10.0));
}

void Meshbuilder::setFixedVertex(int32_t vertex, const Eigen::Vector3d& position)
{
	_V.row(vertex) = position;
}

Eigen::Vector3d Meshbuilder::getFixedVertex(int32_t vertex) const
{
	return _V.row(vertex);
}

void Meshbuilder::gradientTest(double h)
{
	double err = 0.0;
	for (Element* element : _elements)
	{
		err += element->diffTest(_V, h);
	}
	std::cout << "error: " << err << std::endl;
}

void Meshbuilder::computeElements()
{
	// Compute element states
	const int32_t n = _elements.size();
	Eigen::VectorXd energies = Eigen::VectorXd::Zero(n);

	_gradient = Eigen::MatrixXd::Zero(_V.rows(), 3);
	for (int32_t i = 0; i < n; i++)
	{
		Element* element = _elements[i];

		element->compute(_V);
		element->addGradient(_gradient);
		energies[i] = element->getEnergy();
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
			element->colorFaces(lo + alpha * (hi - lo), _C);
		}
	}
	else
	{
		for (int32_t i = 0; i < _elements.size(); i++)
		{
			const Element* element = _elements[i];
			element->colorFaces(lo, _C);
		}
	}
	_energy = energies.sum();
}

void Meshbuilder::gradientDescent(double lambda, int32_t maxIterations)
{
	double step = lambda;

	for (int i = 0; i < maxIterations; i++)
	{
		const auto V = _V;
		const auto G = _gradient;
		const double energy = _energy;
		for (;;)
		{
			_V = V - G * step;
			computeElements();

			if (_energy < energy)
			{
				step *= 1.1;
				break;
			}

			step *= 0.5;

			if (step < EPS)
			{
				_V = V;
				return;
			}
		}
	}
}

void Meshbuilder::renderShell(igl::opengl::glfw::Viewer& viewer) const
{
	viewer.data().clear();
	viewer.data().set_mesh(_V, _F);
	viewer.data().set_colors(_C);
	
	viewer.data().set_face_based(true);
	
	const Eigen::RowVector3d white(1.0, 1.0, 1.0);
	viewer.data().add_edges(_V, _V - _gradient, white);

#ifdef COMPUTE_FINITE_DIFFERENCE
	const Eigen::RowVector3d red(1.0, 0.0, 0.0);
	viewer.data().add_edges(_V, _V - _finiteGradient, red);
#endif
}

/////////////////////////////////////////////////////

void Meshbuilder::computePlane(bool hasCeiling, double altitude)
{
	int32_t n = _oF.rows();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_oV.row(_F(i, j)) = _V.row(_F(i, j));
		}
	}

	Eigen::MatrixXi Ti, Tii;
	igl::triangle_triangle_adjacency(_oF, Ti, Tii);

	//Eigen::MatrixXi Fi, Fii;
	//igl::vertex_triangle_adjacency(_oF, n, Fi, Fii);

	std::vector<double> energies(n);
	std::queue<Face> queue;
	std::set<int32_t> pool;

	// Init all bottom faces and find least deformed
	int root = 0;
	double maxEnergy = 0.0;
	for (int32_t i = 0; i < n; i++)
	{
		pool.insert(i);

		double energy = 0;
		energy += _elements[i * 3 + 0]->getEnergy();
		energy += _elements[i * 3 + 1]->getEnergy();
		energy += _elements[i * 3 + 2]->getEnergy();
		energies[i] = energy;

		if (energy > maxEnergy)
		{
			root = i;
			maxEnergy = energy;
		}
	}

	_pF = Eigen::MatrixXi::Zero(n * (hasCeiling ? 2 : 1), 3);

	// Attach new triangles
	pool.erase(root);
	queue.emplace(Face({ root, 0, -1, -1 }));
	std::vector<Eigen::Vector3d> V;
	while (!queue.empty())
	{
		// Extract next pair
		const Face front = queue.front();
		queue.pop();
		
		// Grab vertices of connected edge
		const auto& face = _oF.row(front.face);
		const Eigen::Vector3i map(
			front.faceEdge,
			(front.faceEdge + 1) % 3,
			(front.faceEdge + 2) % 3);


		Eigen::Vector3d source[3];
		int32_t target[2];
		if (front.prev >= 0)
		{
			// Source connection points
			source[0] = _oV.row(face[map[0]]);
			source[1] = _oV.row(face[map[1]]);
			source[2] = _oV.row(face[map[2]]);
			
			// Target connection points, swap since adjacent faces are differently turned
			const auto prev = _pF.row(front.prev);
			target[1] = prev[front.prevEdge];
			target[0] = prev[(front.prevEdge + 1) % 3];
		}
		else
		{
			source[0] = _oV.row(face[0]);
			source[1] = _oV.row(face[1]);
			source[2] = _oV.row(face[2]);

			// Create anchor points
			target[0] = V.size();
			V.push_back(Eigen::Vector3d(0.0, altitude, 0.0));
			target[1] = V.size();
			V.push_back(Eigen::Vector3d((source[0] - source[1]).norm(), altitude, 0.0));
		}

		// Compute triangle transforms
		Eigen::Matrix3d sourceT;
		const Eigen::Vector3d right = source[2] - source[0];
		sourceT.row(0) = (source[1] - source[0]).normalized();
		sourceT.row(1) = sourceT.row(0).cross(right).normalized();
		sourceT.row(2) = sourceT.row(1).cross(sourceT.row(0));

		Eigen::Matrix3d targetT;
		targetT.row(0) = (V[target[1]] - V[target[0]]).normalized();
		targetT.row(1) = Eigen::Vector3d(0.0, -1.0, 0.0);
		targetT.row(2) = targetT.row(1).cross(targetT.row(0));

		const Eigen::Matrix3d transform = targetT.transpose() * sourceT;

		// Place triangle on plane
		auto& floor = _pF.row(front.face);
		const Eigen::Vector3d ref = V[target[0]];
		floor[map[0]] = target[0];
		floor[map[1]] = target[1];
		floor[map[2]] = V.size();
		V.push_back(ref + (transform * Eigen::Vector3d(source[2] - source[0])));

		if (hasCeiling)
		{
			// Place ceiling
			auto& ceil = _pF.row(front.face + n);
			const auto & refFace = _F.row(front.face + n);
			ceil[map[0]] = V.size();
			V.push_back(ref + (transform * (Eigen::Vector3d(_V.row(refFace[map[0]])) - Eigen::Vector3d(source[0]))));

			ceil[map[1]] = V.size();
			V.push_back(ref + (transform * (Eigen::Vector3d(_V.row(refFace[map[1]])) - Eigen::Vector3d(source[0]))));

			ceil[map[2]] = V.size();
			V.push_back(ref + (transform * (Eigen::Vector3d(_V.row(refFace[map[2]])) - Eigen::Vector3d(source[0]))));
		}

		// Get connected faces
		std::vector<int> raw({0, 1, 2}), edges;

		// Filter for valid faces and order by energy
		std::copy_if(raw.begin(), raw.end(), std::back_inserter(edges), [&](int edge)->bool { return Ti(front.face, edge) >= 0; });
		std::sort(edges.begin(), edges.end(), [&](int a, int b)->bool { return energies[Ti(front.face, a)] < energies[Ti(front.face, b)]; });
		
		for (int edge : edges)
		{
			// only add triangles that aren't already placed
			const int next = Ti(front.face, edge);
			auto it = pool.find(next);
			if (it != pool.end())
			{
				pool.erase(it);
				queue.emplace(Face({ next,  Tii(front.face, edge), front.face, edge }));
			}
		}
	}

	// Create plane
	_pV = Eigen::MatrixXd(V.size(), 3);
	for (int i = 0; i < V.size(); i++)
	{
		_pV.row(i) = V[i];
	}
}

void Meshbuilder::renderPlane(igl::opengl::glfw::Viewer& viewer) const
{
	std::vector<Eigen::MatrixXd> Vs({ _V, _pV });
	std::vector<Eigen::MatrixXi> Fs({ _F, _pF });

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::combine(Vs, Fs, V, F);

	viewer.data().clear();
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
}

/////////////////////////////////////////////////////

IntervalBread Meshbuilder::createIntervalBread(const Eigen::Vector2d& dir, double density) const
{
	// Start with first vertex
	const Eigen::Vector3d& vp = _pV.row(0);
	double min = dir.dot(Eigen::Vector2d(vp.x(), vp.z()));
	for (int i = 0; i < _pF.rows(); i++)
	{
		// Get face corners in slicer space
		const Eigen::Vector3i& face = _pF.row(i);

		const Eigen::Vector3d& va = _pV.row(face.x());
		const double a = dir.dot(Eigen::Vector2d(va.x(), va.z()));

		const Eigen::Vector3d& vb = _pV.row(face.y());
		const double b = dir.dot(Eigen::Vector2d(vb.x(), vb.z()));

		const Eigen::Vector3d& vc = _pV.row(face.z());
		const double c = dir.dot(Eigen::Vector2d(vc.x(), vc.z()));

		min = std::min(a, std::min(b, c));
	}

	// Make slices
	IntervalBread bread(min * dir, dir);
	while (tracePlane(bread, density)) {}
	return bread;
}

bool Meshbuilder::lineIntersection(double y, const Eigen::Vector2d& a, const Eigen::Vector2d& b, Interval& interval) const
{
	if ((a.y() < y) != (b.y() < y))
	{
		const Eigen::Vector2d e = b - a;

		const double d = a.y() - y;
		const double t = d / e.y();
		const double x = a.x() + t * e.x();

		interval.a = std::min(interval.a, x);
		interval.b = std::max(interval.b, x);
		return true;
	}
	return false;
}

bool Meshbuilder::tracePlane(IntervalBread& bread, double density) const
{
	// Make rotation matrix from global to local
	const Eigen::Matrix2d rot = bread.getRot();
	const Eigen::Vector2d pnt = bread.getAnchor();

	IntervalSlice slice;
	const Eigen::Vector2d& p = rot * pnt;
	for (int i = 0; i < _pF.rows(); i++)
	{

		// Get face corners in slicer space
		const Eigen::Vector3i& face = _pF.row(i);

		const Eigen::Vector3d& va = _pV.row(face.x());
		const Eigen::Vector2d& a = rot * Eigen::Vector2d(va.x(), va.z());

		const Eigen::Vector3d& vb = _pV.row(face.y());
		const Eigen::Vector2d& b = rot * Eigen::Vector2d(vb.x(), vb.z());

		const Eigen::Vector3d& vc = _pV.row(face.z());
		const Eigen::Vector2d& c = rot * Eigen::Vector2d(vc.x(), vc.z());

		const double min = std::max(a.x(), std::max(b.x(), c.x()));
		const double max = std::min(a.x(), std::min(b.x(), c.x()));
		Interval interval(min, max);

		bool added = false;
		added = lineIntersection(p.y(), a, b, interval) || added;
		added = lineIntersection(p.y(), a, c, interval) || added;
		added = lineIntersection(p.y(), b, c, interval) || added;
		if (added)
		{
			slice.addInterval(interval);
		}
	}

	// Keep going until nothing has been sliced
	if (!slice.isEmpty())
	{
		bread.walk(density, slice);
		return true;
	}
	return false;
}

void Meshbuilder::renderSlices(const IntervalBread& bread, igl::opengl::glfw::Viewer& viewer, double altitude) const
{
	std::vector<Eigen::Vector2d> Lfrom;
	std::vector<Eigen::Vector2d> Lto;
	bread.depthSearch([&](const Eigen::Vector2d & from, const Eigen::Vector2d & to)
	{
			Lfrom.push_back(from);
			Lto.push_back(to);
	});

	Eigen::MatrixXd VFrom(Lfrom.size(), 3);
	Eigen::MatrixXd VTo(Lto.size(), 3);
	for (int i = 0; i < VFrom.rows(); i++)
	{
		VFrom.row(i) = Eigen::Vector3d(Lfrom[i].x(), altitude, Lfrom[i].y());
		VTo.row(i) = Eigen::Vector3d(Lto[i].x(), altitude, Lto[i].y());
	}

	const Eigen::RowVector3d red(1.0, 0.0, 0.0);
	viewer.data().add_edges(VFrom, VTo, red);
}