
#pragma once

#include <Eigen/Core>
#include <vector>
#include <functional>

#define INTERVAL_EPS 0.001


struct Interval
{
	Interval(double a, double b);

	bool operator<(double x) const;
	bool operator<(const Interval& other) const;
	bool operator>(double x) const;
	bool operator>(const Interval& other) const;
	bool inside(double x) const;
	bool inside(const Interval& other) const;

	bool combine(const Interval& other);

	double a;
	double b;
};


class IntervalSlice
{
	friend class IntervalBread;

	/////////////////////////////////////////////////////
public:

	IntervalSlice();
	~IntervalSlice();

	void addInterval(const Interval& interval);
	bool isEmpty() const;

private:

	std::vector<Interval> _intervals;

};


class IntervalBread
{
	/////////////////////////////////////////////////////
public:

	IntervalBread(const Eigen::Vector2d& org, const Eigen::Vector2d& dir);
	~IntervalBread();

	void walk(double distance, const IntervalSlice& slice);
	Eigen::Vector2d getAnchor() const;
	Eigen::Matrix2d getRot() const;

	void depthSearch(std::function<void(const Eigen::Vector2d& from, const Eigen::Vector2d& to)> func) const;

private:
	Eigen::Matrix2d _rot;
	Eigen::Vector2d _org;
	double _distance;
	std::vector<IntervalSlice> _slices;

};