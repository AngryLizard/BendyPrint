#include "meshunwrapper.h"
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/combine.h>

#include <set>
#include <queue>

MeshUnwrapper::MeshUnwrapper(const Eigen::MatrixXd& oV, const Eigen::MatrixXi& oF, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double altitude, std::function<bool(int32_t, int32_t)> canConnect)
{
	int32_t n = oF.rows();

	// Adjacency is used to know which faces to place next to each other
	Eigen::MatrixXi Ti, Tii;
	igl::triangle_triangle_adjacency(oF, Ti, Tii);
	
	// All faces are pooled so none is forgotten
	std::set<int32_t> pool;

	// TODO: Implement a smart way of determining root face
	int root = 0;
	for (int i = 0; i < oF.rows(); i++)
	{
		pool.emplace(i);
	}

	// Init plane vertices/faces
	_pF = Eigen::MatrixXi::Zero(n * 1, 3);
	_dF = Eigen::MatrixXi::Zero(n * 1, 3);

	// Attach new triangles
	pool.erase(root);

	// Structure to keep track of queued face attachements
	struct Face
	{
		int32_t face; // Face to be projected
		int32_t faceEdge; // Face edge to attach
		int32_t prev; // Face to be attached to
		int32_t prevEdge; // Edge to be attach to
	};
	std::queue<Face> queue;
	queue.emplace(Face({ root, 0, -1, -1 }));

	// New faces/vertices are queued up in vectors since final number isn't known yet for vertices
	std::vector<Eigen::Vector3d> pV, dV;
	while (!queue.empty())
	{
		// Extract next face to attach from queue
		const Face front = queue.front();
		queue.pop();

		// Grab vertices of connected edge, mapping determines which vertices are touching
		const auto& face = oF.row(front.face);
		const Eigen::Vector3i map(
			front.faceEdge,
			(front.faceEdge + 1) % 3,
			(front.faceEdge + 2) % 3);

		// Three -source- vertices are attached to two -target- vertices depending on which vertices are touching on original mesh
		Eigen::Vector3d source[3];
		int32_t target[2];
		if (front.prev >= 0)
		{
			// Source connection points
			source[0] = oV.row(face[map[0]]);
			source[1] = oV.row(face[map[1]]);
			source[2] = oV.row(face[map[2]]);

			// Target connection points
			const auto prev = _pF.row(front.prev);
			target[1] = prev[front.prevEdge];
			target[0] = prev[(front.prevEdge + 1) % 3];
		}
		else
		{
			// Initial face gets put on 0,0
			source[0] = oV.row(face[0]);
			source[1] = oV.row(face[1]);
			source[2] = oV.row(face[2]);

			// Create anchor points
			target[0] = pV.size();
			pV.push_back(Eigen::Vector3d(0.0, altitude, 0.0));
			target[1] = pV.size();
			pV.push_back(Eigen::Vector3d((source[0] - source[1]).norm(), altitude, 0.0));
		}

		// Compute triangle transforms to transforms other vertices
		Eigen::Matrix3d sourceT;
		const Eigen::Vector3d right = source[2] - source[0];
		sourceT.row(0) = (source[1] - source[0]).normalized();
		sourceT.row(1) = sourceT.row(0).cross(right).normalized();
		sourceT.row(2) = sourceT.row(1).cross(sourceT.row(0));

		Eigen::Matrix3d targetT;
		targetT.row(0) = (pV[target[1]] - pV[target[0]]).normalized();
		targetT.row(1) = Eigen::Vector3d(0.0, -1.0, 0.0);
		targetT.row(2) = targetT.row(1).cross(targetT.row(0));

		const Eigen::Matrix3d transform = targetT.transpose() * sourceT;

		// Place triangle on plane so the third vertex matches
		auto& floor = _pF.row(front.face);
		const Eigen::Vector3d ref = pV[target[0]];
		floor[map[0]] = target[0];
		floor[map[1]] = target[1];
		floor[map[2]] = pV.size();
		pV.push_back(ref + (transform * Eigen::Vector3d(source[2] - source[0])));

		// Also place ceiling according to previously computed transform
		auto& ceil = _dF.row(front.face);
		const auto& refFace = F.row(front.face + n);

		ceil[map[0]] = dV.size();
		dV.push_back(ref + (transform * (Eigen::Vector3d(V.row(refFace[map[0]])) - Eigen::Vector3d(source[0]))));

		ceil[map[1]] = dV.size();
		dV.push_back(ref + (transform * (Eigen::Vector3d(V.row(refFace[map[1]])) - Eigen::Vector3d(source[0]))));

		ceil[map[2]] = dV.size();
		dV.push_back(ref + (transform * (Eigen::Vector3d(V.row(refFace[map[2]])) - Eigen::Vector3d(source[0]))));

		// Get connected faces to newly added face to add to the queue
		std::vector<int> raw({ 0, 1, 2 }), edges;

		// Filter for valid face connections
		std::copy_if(raw.begin(), raw.end(), std::back_inserter(edges), [&](int edge)->bool { return Ti(front.face, edge) >= 0 && canConnect(front.face, Ti(front.face, edge)); });
		for (int edge : edges)
		{
			// Only add triangles to the queue that aren't already placed
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
	_pV = Eigen::MatrixXd(pV.size(), 3);
	for (int i = 0; i < pV.size(); i++)
	{
		_pV.row(i) = pV[i];
	}

	const double mergeThreshold = 0.1;
	auto replaceVertex = [&](int32_t from, int32_t to) {

		for (int i = 0; i < _pF.rows(); i++)
		{
			auto& face = _pF.row(i);
			if (face.x() == from) face.x() = to;
			if (face.y() == from) face.y() = to;
			if (face.z() == from) face.z() = to;
		}
	};

	// Remove duplicates (merges connected edges)
	for (int i = 0; i < _pF.rows(); i++)
	{
		const auto& face = _pF.row(i);
		for (int j = i + 1; j < _pF.rows(); j++)
		{
			auto& other = _pF.row(j);
			if ((_pV.row(face.x()) - _pV.row(other.x())).norm() < mergeThreshold) replaceVertex(other.x(), face.x());
			if ((_pV.row(face.y()) - _pV.row(other.y())).norm() < mergeThreshold) replaceVertex(other.y(), face.y());
			if ((_pV.row(face.z()) - _pV.row(other.z())).norm() < mergeThreshold) replaceVertex(other.z(), face.z());
		}
	}

	// Create ceiling
	_dV = Eigen::MatrixXd(dV.size(), 3);
	for (int i = 0; i < dV.size(); i++)
	{
		_dV.row(i) = dV[i];
	}

	_slicer = new Slicer(_pV, _pF, _dV, _dF);
}

MeshUnwrapper::~MeshUnwrapper()
{
	if (_slicer)
	{
		delete(_slicer);
	}
}

void MeshUnwrapper::renderPlane(igl::opengl::glfw::Viewer& viewer, bool hasCeiling) const
{
	std::vector<Eigen::MatrixXd> Vs({ _pV });
	std::vector<Eigen::MatrixXi> Fs({ _pF });

	if (hasCeiling)
	{
		Vs.push_back(_dV);
		Fs.push_back(_dF);
	}

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::combine(Vs, Fs, V, F);

	viewer.data().clear();
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
}

/////////////////////////////////////////////////////

const Slicer* MeshUnwrapper::getSlicer() const
{
	return _slicer;
}
