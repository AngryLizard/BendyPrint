
#include "slicer.h"
#include <iostream>

#include <igl/triangle_triangle_adjacency.h>


Slicer::Slicer(const Eigen::MatrixXd& planeV, const Eigen::MatrixXi& planeF, const Eigen::MatrixXd& fillV, const Eigen::MatrixXi& fillF)
	: _planeV(planeV), _planeF(planeF), _fillV(fillV), _fillF(fillF)
{
}

Slicer::~Slicer()
{
}


Eigen::Vector2d Slicer::flatten(const Eigen::Vector3d& vec) const
{
	return Eigen::Vector2d(vec.x(), vec.z());
}

Eigen::Vector3d Slicer::project(const Eigen::Vector2d& vec, double altitude) const
{
	return Eigen::Vector3d(vec.x(), altitude, vec.y());
}



IntervalBread Slicer::createIntervalBread(const Eigen::Vector2d& dir, double margin) const
{
	// Start with first vertex
	double min = FLT_MAX;
	for (int i = 0; i < _planeV.rows(); i++)
	{
		// Get face corners in slicer space
		const Eigen::Vector3d& va = _planeV.row(i);
		const double a = dir.dot(Eigen::Vector2d(va.x(), va.z()));
		min = std::min(min, a);
	}

	// Make slices
	IntervalBread bread(min * dir, dir);
	while (tracePlane(bread, margin)) {}
	return bread;
}

bool Slicer::lineIntersection(const Eigen::Vector2d& a, const Eigen::Vector2d& b, Interval& interval) const
{
	if (((a.x() < EPS) && (b.x() > -EPS)) || ((a.x() > -EPS) && (b.x() < EPS)))
	{
		const Eigen::Vector2d e = a - b;
		if (e.x() * e.x() > EPS)
		{
			const double t = a.x() / e.x();
			const double y = a.y() - t * e.y();
			interval.a = std::min(interval.a, y);
			interval.b = std::max(interval.b, y);
		}
		else
		{
			interval.a = std::min(interval.a, std::min(a.y(), b.y()));
			interval.b = std::max(interval.b, std::max(a.y(), b.y()));
		}

		return true;
	}
	return false;
}

bool Slicer::tracePlane(IntervalBread& bread, double margin) const
{
	// Make rotation matrix from global to local
	const Eigen::Matrix2d rot = bread.getRot();
	const Eigen::Vector2d pnt = bread.getCurrent();

	IntervalSlice slice;
	for (int i = 0; i < _planeF.rows(); i++)
	{

		// Get face corners in slicer space
		const Eigen::Vector3i& face = _planeF.row(i);

		const Eigen::Vector2d& a = rot * (flatten(_planeV.row(face.x())) - pnt);
		const Eigen::Vector2d& b = rot * (flatten(_planeV.row(face.y())) - pnt);
		const Eigen::Vector2d& c = rot * (flatten(_planeV.row(face.z())) - pnt);

		const double min = std::max(a.y(), std::max(b.y(), c.y()));
		const double max = std::min(a.y(), std::min(b.y(), c.y()));
		Interval interval(min, max);

		bool added = false;
		added = lineIntersection(a, b, interval) || added;
		added = lineIntersection(a, c, interval) || added;
		added = lineIntersection(b, c, interval) || added;
		if (added)
		{
			slice.addInterval(interval);
		}
	}

	// Keep going until nothing has been sliced
	if (!slice.isEmpty())
	{
		bread.walk(margin, slice);
		return true;
	}
	return false;
}


IntervalBread Slicer::createIntervalBread(const Eigen::Vector2d& dir, const OutlineBread& outlines, double margin) const
{
	// Find first vertex location
	double min = FLT_MAX;
	outlineSearch(outlines, [&min, dir](const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool)
		{
			const double x = std::min(dir.dot(from), dir.dot(to));
			min = std::min(min, x);
		}, margin);

	// Make slices
	IntervalBread bread(min * dir, dir);
	while (tracePlane(bread, outlines, margin)) {}
	return bread;
}

bool Slicer::tracePlane(IntervalBread& bread, const OutlineBread& outlines, double margin) const
{
	using LocalPair = std::pair< Eigen::Vector2d, Eigen::Vector2d>;
	using Intersection = std::pair<double, bool>; // Y coordinate / X sign bit

	// Make rotation matrix from global to local
	const Eigen::Matrix2d rot = bread.getRot();
	const Eigen::Vector2d pnt = bread.getCurrent();

	std::vector<LocalPair> localLines;
	outlineSearch(outlines, [&localLines, rot, pnt](const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool)
		{
			localLines.push_back(LocalPair(rot * (from - pnt), rot * (to - pnt)));
		}, margin);

	// Get all intersections
	std::vector<Intersection> intersections;
	for (const LocalPair& pair : localLines)
	{
		double point;
		if (((pair.first.x() < EPS) && (pair.second.x() > -EPS)) || ((pair.first.x() > -EPS) && (pair.second.x() < EPS)))
		{
			const Eigen::Vector2d e = pair.first - pair.second;
			if (e.x() * e.x() > EPS)
			{
				const double t = pair.first.x() / e.x();
				point = pair.first.y() - t * e.y();
				
				// Keep track of the sign of the further end to detect intersections on edge cases
				const int32_t sign = (t < 0.5) ? (pair.second.x() < 0.0) : (pair.first.x() < 0.0);
				intersections.push_back(Intersection(point, sign));
			}
			else
			{
				intersections.push_back(Intersection(pair.first.x(), pair.second.x() < 0.0));
				intersections.push_back(Intersection(pair.second.x(), pair.first.x() < 0.0));
			}
		}
	}

	// By sorting, every two points should make an interval inside the mesh.
	std::sort(intersections.begin(), intersections.end(), [](const Intersection& a, const Intersection& b) { return a.first < b.first; });

	// Remove duplicates
	if (intersections.size() > 2)
	{
		auto it = intersections.begin();
		while (it + 1 < intersections.end())
		{
			const Intersection its = *(it + 1);
			const double dist = its.first - it->first;
			if (dist * dist < EPS)
			{
				// On a corner we can remove the intersection altogether
				if (its.second == it->second)
				{
					it = intersections.erase(it);
				}
				it = intersections.erase(it);
			}
			else
			{
				it++;
			}
		}
	}

	// Keep going until nothing has been sliced
	if (!intersections.empty())
	{
		assert((intersections.size()&1) == 0 && "Outline trace always returns an even amount of intersections.");

		IntervalSlice slice;
		for (int i = 0; i < intersections.size(); i += 2)
		{
			slice.addInterval(Interval(intersections[i].first, intersections[i + 1].first));
		}
		bread.walk(margin, slice);
		return true;
	}
	return false;
}


OutlineBread Slicer::createOutlineBread() const
{
	Eigen::MatrixXi Ti, Tii;
	igl::triangle_triangle_adjacency(_planeF, Ti, Tii);

	OutlineBread bread;
	for (int i = 0; i < _planeF.rows(); i++)
	{
		const Eigen::Vector3i& face = _planeF.row(i);

		// Add edges without neighbour
		int32_t neighbourCount = 0;
		for (int j = 0; j < 3; j++)
		{
			if (Ti(i, j) == -1)
			{
				bread.addOutline(OutlineSlice(face[j], face[(j + 1) % 3]));
			}
		}
	}

	return bread;
}

void Slicer::outlineSearch(const OutlineBread& bread, std::function<void(const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)> func, double margin) const
{
	bool jump = false;
	int32_t currentIndex = -1;
	Eigen::Vector2d first, last;
	bread.depthSearch([&](const Eigen::Vector3i& line, int32_t index)
		{
			// Compute offset point with given margin
			const Eigen::Vector2d prev = flatten(_planeV.row(line.x()));
			const Eigen::Vector2d curr = flatten(_planeV.row(line.y()));
			const Eigen::Vector2d next = flatten(_planeV.row(line.z()));

			const Eigen::Vector2d diff = curr - prev;
			const Eigen::Vector2d diffOrth = Eigen::Vector2d(diff.y(), -diff.x()).normalized();

			const Eigen::Vector2d gap = (next - curr).normalized() + diff.normalized();
			const Eigen::Vector2d gapOrth = Eigen::Vector2d(gap.y(), -gap.x()).normalized();

			// Inv should always be positive if no edges overlap
			const double inv = diffOrth.dot(gapOrth);
			if (inv > EPS)
			{
				const bool swtch = (currentIndex != index);
				if (currentIndex > -1 && swtch)
				{
					// Close loop
					func(last, first, false);
				}

				// Make curve if too pointy
				const Eigen::Vector2d gapNormal = gap.normalized();
				const double dot = diffOrth.dot(gapNormal);
				if (dot > EPS)
				{
					const Eigen::Vector2d target = curr - gapOrth * margin;
					const Eigen::Vector2d offset = gapNormal * margin * dot;
					const Eigen::Vector2d targetA = target - offset;
					const Eigen::Vector2d targetB = target + offset;

					if (swtch)
					{
						currentIndex = index;
						first = targetA;
					}
					else
					{
						func(last, targetA, jump);
						jump = false;
					}

					func(targetA, targetB, swtch);
					last = targetB;

				}
				else
				{
					// Offset orthogonally and remember for next point
					const Eigen::Vector2d target = curr - gapOrth * (margin / inv);
					if (swtch)
					{
						currentIndex = index;
						first = target;
						jump = true;
					}
					else
					{
						func(last, target, jump);
						jump = false;
					}

					last = target;
				}
			}
		});

	if (currentIndex > -1)
	{
		// Close loop
		func(last, first, false);
	}
}



TriangleBread Slicer::createTriangleBread(double scale, double height) const
{
	TriangleBread bread;
	Eigen::MatrixXi Ti, Tii;
	igl::triangle_triangle_adjacency(_planeF, Ti, Tii);

	// Add all faces to triangle pool
	const int32_t faces = std::min(_planeF.rows(), _fillF.rows());
	for (int i = 0; i < faces; i++)
	{
		const Eigen::Vector3i& floor = _planeF.row(i);
		const Eigen::Vector3i& ceil = _fillF.row(i);

		double lower = 0.0;
		double upper = 0.0;

		// Go through all vertices and store final positions
		Eigen::Vector2d floors[3];
		Eigen::Vector2d ceilings[3];
		for (int j = 0; j < 3; j++)
		{
			const int32_t a = (j + 1) % 3;
			const int32_t b = (j + 2) % 3;

			// Position primary vertex according to whether neighbouring edges are free or connected
			floors[j] = flatten(_planeV.row(floor[j]));
			if (Ti(i, j) != -1 || Ti(i, b) != -1)
			{
				const Eigen::Vector2d bSecondary = flatten(_planeV.row(floor[a]));
				const Eigen::Vector2d bTernary = flatten(_planeV.row(floor[b]));

				// Put a gap on connected edges
				const Eigen::Vector2d bEdge = bTernary - bSecondary;
				const double alpha = (Ti(i, j) == -1 ? 0.0 : (Ti(i, b) == -1 ? 1.0 : 0.5));
				const Eigen::Vector2d bCenter = (bSecondary + bEdge * alpha);
				floors[j] = bCenter + (floors[j] - bCenter) * (1.0 + scale) / 2;
			}
			ceilings[j] = flatten(_fillV.row(ceil[j]));

			lower += _planeV.row(floor[j]).y();
			upper += _fillV.row(ceil[j]).y();
		}

		// Swap x and y to make triangulation coherent
		std::swap(floors[0], floors[1]);

		// Interpolate vertices according to depth
		const double depth = height / ((upper - lower) / 3);
		Eigen::Vector2d vertices[3];
		for (int j = 0; j < 3; j++)
		{
			vertices[j] = floors[j] + depth * (ceilings[j] - floors[j]);
		}

		bread.addTriangle(vertices);
	}

	bread.sortTriangles();
	return bread;
}

void Slicer::triangleFillSearch(const TriangleBread& bread, std::function<void(const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)> func, int32_t layers, double margin) const
{
	bread.depthSearch([&](const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool fresh) {
		func(from, to, fresh);
		}, layers, margin);
}
