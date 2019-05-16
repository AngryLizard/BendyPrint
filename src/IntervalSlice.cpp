
#include "IntervalSlice.h"

Interval::Interval(double a, double b)
: a(std::min(a, b)), b(std::max(a, b))
{
}

bool Interval::operator<(double x) const
{
	return x <= a;
}

bool Interval::operator>(double x) const
{
	return x >= b;
}

bool Interval::inside(double x) const
{
	return a <= x && x <= b;
}

bool Interval::operator<(const Interval& other) const
{
	return b <= other.a;
}

bool Interval::operator>(const Interval& other) const
{
	return a >= other.b;
}

bool Interval::inside(const Interval& other) const
{
	return other.a <= a && b <= other.b;
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


IntervalSlice::IntervalSlice()
{

}

IntervalSlice::~IntervalSlice()
{

}


void IntervalSlice::addInterval(const Interval& interval)
{

	auto it = _intervals.begin();

	// Go until overlap
	int32_t idx = 0;
	while (it != _intervals.end() && interval < *it)
	{
		it++;
	}

	if (it == _intervals.end())
	{
		_intervals.push_back(interval);
	}
	else
	{
		// Integrate if possible, insert otherwise
		if (!it->combine(interval))
		{
			_intervals.insert(it, interval);
		}

		// Check whether future intervals can be combined
		while ((it + 1) != _intervals.end())
		{
			it++;
			if (it->combine(*(it - 1)))
			{
				_intervals.erase(it);
			}
		}
	}
}

bool IntervalSlice::isEmpty() const
{
	return _intervals.empty();
}


IntervalBread::IntervalBread(const Eigen::Vector2d& org, const Eigen::Vector2d& dir)
	: _org(org), _distance(0.0)
{
	_rot.row(0) = Eigen::Vector2d(-dir.y(), dir.x());
	_rot.row(1) = dir;
}

IntervalBread::~IntervalBread()
{

}

void IntervalBread::walk(double distance, const IntervalSlice& slice)
{
	_slices.push_back(slice);
	_distance += distance;
	_org += _rot * Eigen::Vector2d(0.0, _distance);
}

Eigen::Vector2d IntervalBread::getAnchor() const
{
	return _org;
}

Eigen::Matrix2d IntervalBread::getRot() const
{
	return _rot;
}

void IntervalBread::depthSearch(std::function<void(const Eigen::Vector2d&, const Eigen::Vector2d&)> func) const
{
	// Depth search through all slices
	int32_t branch = 0;
	bool found;
	do
	{
		// Go until full depth was visited
		found = false;
		for (int i = 0; i < _slices.size(); i++)
		{
			const IntervalSlice& slice = _slices[i];
			if (branch < slice._intervals.size())
			{
				const Interval& interval = slice._intervals[branch];
				const double y = (_distance / _slices.size()) * i;
				const Eigen::Vector2d from = _org + _rot * Eigen::Vector2d(interval.a, y);
				const Eigen::Vector2d to = _org + _rot * Eigen::Vector2d(interval.b, y);

				func(from, to);
				found = true;
			}
		}
		branch++;
	} while (found);
}