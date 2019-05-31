
#include "bread.h"
#include <iostream>
#include <queue>

/////////////////////////////////////////////////////

Interval::Interval(double a, double b)
: a(a), b(b)
{
}

bool Interval::operator<(double x) const
{
	return x <= a - INTERVAL_EPS;
}

bool Interval::operator>(double x) const
{
	return x >= b + INTERVAL_EPS;
}

bool Interval::inside(double x) const
{
	return a <= x - INTERVAL_EPS && x + INTERVAL_EPS <= b;
}

bool Interval::operator<(const Interval& other) const
{
	return b <= other.a - INTERVAL_EPS;
}

bool Interval::operator>(const Interval& other) const
{
	return a >= other.b + INTERVAL_EPS;
}

bool Interval::inside(const Interval& other) const
{
	return other.a <= a - INTERVAL_EPS && b + INTERVAL_EPS <= other.b;
}


bool Interval::combine(const Interval& other)
{
	if (b <= other.a - INTERVAL_EPS || a >= other.b + INTERVAL_EPS)
	{
		return false;
	}

	a = std::min(a, other.a);
	b = std::max(b, other.b);

	return true;
}

double Interval::area(const Interval& other) const
{
	if (other < *this || *this < other)
	{
		return 0.0;
	}


	const double ta = std::max(a, other.a);
	const double tb = std::min(b, other.b);
	return tb - ta;
}

/////////////////////////////////////////////////////

IntervalSlice::IntervalSlice()
{

}

IntervalSlice::~IntervalSlice()
{

}


void IntervalSlice::addInterval(const Interval& interval)
{
	// See whether interval can be combined with existing interval
	for (auto it = _intervals.begin(); it < _intervals.end(); it++)
	{
		if (it->combine(interval))
		{
			// Remove from list and add it back until no combinations are possible anymore
			const Interval combined = *it;
			it = _intervals.erase(it);
			addInterval(combined);
			return;
		}
	}

	// Add interval and spatially sort
	_intervals.push_back(interval);
	std::sort(_intervals.begin(), _intervals.end(), [](const Interval & a, const Interval & b)->bool { return a < b; });
}

bool IntervalSlice::isEmpty() const
{
	return _intervals.empty();
}

/////////////////////////////////////////////////////

IntervalBread::IntervalBread(const Eigen::Vector2d& org, const Eigen::Vector2d& dir)
	: _org(org), _distance(0.0)
{
	const Eigen::Vector2d& nrm = dir.normalized();
	_rot.row(0) = nrm;
	_rot.row(1) = Eigen::Vector2d(-nrm.y(), nrm.x());
}

IntervalBread::~IntervalBread()
{

}

void IntervalBread::walk(double distance, const IntervalSlice& slice)
{
	_slices.push_back(slice);
	_distance += distance;
}

Eigen::Vector2d IntervalBread::getCurrent() const
{
	return _org + _rot.transpose() * Eigen::Vector2d(_distance, 0.0);
}

Eigen::Vector2d IntervalBread::getAnchor() const
{
	return _org;
}

Eigen::Matrix2d IntervalBread::getRot() const
{
	return _rot;
}


void IntervalBread::depthSearch(std::function<void(const Eigen::Vector2d&, const Eigen::Vector2d&, bool)> func, double margin, double threshold) const
{
	using Pair = std::pair<Eigen::Vector2d, Eigen::Vector2d>;

	// store sequence for reordering
	std::vector<std::vector<Pair>> sequence;

	// Make working copy and go through all slices until all intervals are accounted for
	std::vector<IntervalSlice> slices = _slices;
	for (int i = 0; i < slices.size(); i++)
	{
		IntervalSlice& base = slices[i];
		while (!base._intervals.empty())
		{
			// Pop first interval
			auto it = base._intervals.begin();
			Interval interval = *it;
			it = base._intervals.erase(it);

			// Transform first interval to 2D space for processing, with jumping
			const double x = (_distance / _slices.size()) * i;
			if (interval.b - interval.a > 2 * margin)
			{
				const Eigen::Vector2d from = _org + _rot.transpose() * Eigen::Vector2d(x, interval.a + margin);
				const Eigen::Vector2d to = _org + _rot.transpose() * Eigen::Vector2d(x, interval.b - margin);
				sequence.push_back(std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>(1, Pair(from, to)));
				//func(from, to, true);
			}
			else
			{
				// Even if slice is too small, always process first slice
				const Eigen::Vector2d from = _org + _rot.transpose() * Eigen::Vector2d(x, interval.a);
				const Eigen::Vector2d to = _org + _rot.transpose() * Eigen::Vector2d(x, interval.b);
				const Eigen::Vector2d avg = (from + to) / 2;
				sequence.push_back(std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>(1, Pair(avg, avg)));
				//func(avg, avg, true);
			}

			// Go through future slices to combine common areas
			for (int j = i + 1; j < slices.size(); j++)
			{
				IntervalSlice& slice = slices[j];

				// Determine area for thresholding
				const double area = interval.b - interval.a;

				// Look for matching interval
				bool found = false;
				for (auto pt = slice._intervals.begin(); pt != slice._intervals.end(); pt++)
				{
					const double match = pt->area(interval) / area;
					const bool areaMatch = (match > threshold && match < 1.0 / threshold);
					if (pt->inside(interval) || interval.inside(*pt) || areaMatch)
					{
						// Transform interval to 2D space for processing, without jumping
						if (pt->b - pt->a > 2 * margin)
						{
							const double x = (_distance / _slices.size()) * j;
							const Eigen::Vector2d from = _org + _rot.transpose() * Eigen::Vector2d(x, pt->a + margin);
							const Eigen::Vector2d to = _org + _rot.transpose() * Eigen::Vector2d(x, pt->b - margin);
							Pair(from, to);

							// Find closer end from last sequence
							const Pair& last = sequence.back().back();
							if ((last.second - from).squaredNorm() < (last.second - to).squaredNorm())
							{
								sequence.back().push_back(Pair(from, to));
							}
							else
							{
								sequence.back().push_back(Pair(to, from));
							}
							//func(from, to, false);
						}

						// Pop interval, remember it for the next slice und go to next slice
						interval = *pt;
						pt = slice._intervals.erase(pt);
						found = true;
						break;
					}
				}

				// Start new sequence if no matching slice was found
				if (!found)
				{
					break;
				}
			}
		}
	}

	/*
	TODO: Reorder to minimize traveling distance
	struct { int32_t a, b; double distance; } best;
	best.distance = (sequence.front().back().first - sequence.back().front().first).squaredNorm();
	for (int i = 0; i < sequence.size(); i++)
	{
		for (int j = 0; j < sequence.size(); j++)
		{
			if (i != j)
			{
				const double a = (sequence[i].back().first - sequence[j].front().first).squaredNorm();
				const double b = (sequence[i].back().first - sequence[j].front().second).squaredNorm();
				const double c = (sequence[i].back().second - sequence[j].front().first).squaredNorm();
				const double d = (sequence[i].back().second - sequence[j].front().second).squaredNorm();
				const double distance = std::min(std::min(a, b), std::min(c, d));
				if (distance < best.distance)
				{
					best.distance = distance;
					best.a = i;
					best.b = j;
				}
			}
		}
	}
	*/
	
	// Play back sequence
	for (const std::vector<Pair>& sub : sequence)
	{
		func(sub.front().first, sub.front().second, true);
		for (int i = 1; i < sub.size(); i++)
		{
			func(sub[i].first, sub[i].second, false);
		}
	}
}

/////////////////////////////////////////////////////

OutlineSlice::OutlineSlice(int32_t a, int32_t b)
{
	_vertices.resize(2, -1);
	_vertices[0] = a;
	_vertices[1] = b;
}

OutlineSlice::~OutlineSlice()
{
}

bool OutlineSlice::combine(const OutlineSlice& other)
{
	if (_vertices.front() == other._vertices.back())
	{
		for (int i = other._vertices.size() - 2; i >= 0; i--)
		{
			_vertices.push_front(other._vertices[i]);
		}
		return true;
	}
	else if (_vertices.front() == other._vertices.front())
	{
		for (int i = 1; i < other._vertices.size(); i++)
		{
			_vertices.push_front(other._vertices[i]);
		}
		return true;
	}
	else if (_vertices.back() == other._vertices.front())
	{
		for (int i = 1; i < other._vertices.size(); i++)
		{
			_vertices.push_back(other._vertices[i]);
		}
		return true;
	}
	else if (_vertices.back() == other._vertices.back())
	{
		for (int i = other._vertices.size() - 2; i >= 0; i--)
		{
			_vertices.push_back(other._vertices[i]);
		}
		return true;
	}
	return false;
}

/////////////////////////////////////////////////////

OutlineBread::OutlineBread()
{
}

OutlineBread::~OutlineBread()
{
}

void OutlineBread::addOutline(const OutlineSlice& outline)
{
	for (auto it = _slices.begin(); it < _slices.end(); it++)
	{
		if (it->combine(outline))
		{
			// Remove from list and add it back until no combinations are possible anymore
			const OutlineSlice combined = *it;
			it = _slices.erase(it);
			addOutline(combined);
			return;
		}
	}
	_slices.push_back(outline);
}


void OutlineBread::depthSearch(std::function<void(const Eigen::Vector3i& line, int32_t index)> func) const
{
	for (int i = 0; i < _slices.size(); i++)
	{
		const OutlineSlice& slice = _slices[i];

		const int32_t n = slice._vertices.size();
		if (n > 1)
		{
			// Last element is the same as the first
			func(Eigen::Vector3i(slice._vertices[n - 2], slice._vertices[0], slice._vertices[1]), i);
			for (int j = 1; j < n-1; j++)
			{
				func(Eigen::Vector3i(slice._vertices[j - 1], slice._vertices[j], slice._vertices[j + 1]), i);
			}
		}
	}
}

/////////////////////////////////////////////////////

Triangle::Triangle(const Eigen::Vector2d vertices[3])
{
	_vertices[0] = vertices[0];
	_vertices[1] = vertices[1];
	_vertices[2] = vertices[2];

	_center = (_vertices[0] + _vertices[1] + _vertices[2]) / 3;
}

Triangle::~Triangle()
{
}



TriangleBread::TriangleBread()
{
}

TriangleBread::~TriangleBread()
{
}

void TriangleBread::addTriangle(const Eigen::Vector2d vertices[3])
{
	_triangles.push_back(Triangle(vertices));
}

void TriangleBread::sortTriangles()
{
	std::vector<Triangle> triangles(_triangles);
	_triangles.clear();
	_triangles.push_back(triangles.back());
	triangles.pop_back();

	while (!triangles.empty())
	{
		// Take next triangle
		_triangles.push_back(triangles.back());
		triangles.pop_back();

		// Sort by proximity (TODO: Replace by priority queue)
		std::sort(triangles.begin(), triangles.end(), [&](const Triangle& a, const Triangle& b)->bool { return (a._center - _triangles.back()._center).squaredNorm() > (b._center - _triangles.back()._center).squaredNorm(); });
	}
}

void TriangleBread::depthSearch(std::function<void(const Eigen::Vector2d&, const Eigen::Vector2d&, bool)> func, int32_t layers, double margin) const
{
	for (const Triangle& triangle : _triangles)
	{
		for (int i = 0; i < layers; i++)
		{
			// Create triangles with given margins
			Eigen::Vector2d vertices[3];
			for (int j = 0; j < 3; j++)
			{
				const Eigen::Vector2d primary = triangle._vertices[j];
				const Eigen::Vector2d secondary = triangle._vertices[(j + 1) % 3];
				const Eigen::Vector2d ternary = triangle._vertices[(j + 2) % 3];

				const Eigen::Vector2d aEdge = (secondary - primary);
				const double aNorm = aEdge.norm();
				const Eigen::Vector2d aNormal = aEdge / aNorm;
				const Eigen::Vector2d aOrth = Eigen::Vector2d(aNormal.y(), -aNormal.x());

				const Eigen::Vector2d bEdge = (ternary - primary);
				const double bNorm = bEdge.norm();
				const Eigen::Vector2d bNormal = bEdge / bNorm;
				const Eigen::Vector2d bOrth = Eigen::Vector2d(-bNormal.y(), bNormal.x());

				const double marg = margin * ((double)i + 0.5);
				const Eigen::Vector2d margins(marg / aOrth.dot(bNormal), marg / bOrth.dot(aNormal));
				
				if (margins.x() / bNorm > 0.5 || margins.y() / aNorm > 0.5)
				{
					// Abort if triangle is filled
					break;
				}

				vertices[j] = primary + bNormal * margins.x() + aNormal * margins.y();
			}

			// Render out all edges
			func(vertices[0], vertices[1], i == 0);
			func(vertices[1], vertices[2], false);
			func(vertices[2], vertices[0], false);
		}
	}
}