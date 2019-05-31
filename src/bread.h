
#pragma once

///////////////////////////////////////////////////////////////
// I call them bread because they have slices don't judge me //
///////////////////////////////////////////////////////////////

#include <Eigen/Core>

#include <set>
#include <deque>
#include <vector>
#include <functional>

#define INTERVAL_EPS 0.001

// Interval structure with common operations
struct Interval
{
	Interval(double a, double b);

	bool operator<(double x) const;
	bool operator<(const Interval& other) const;
	bool operator>(double x) const;
	bool operator>(const Interval& other) const;
	bool inside(double x) const;
	bool inside(const Interval& other) const;

	// Combines if possible, return true if it did
	bool combine(const Interval& other);

	// Area between two intervals, 0 if not intersecting
	double area(const Interval& other) const;

	double a;
	double b;
};

// Collection of intervals
class IntervalSlice
{
	friend class IntervalBread;

	/////////////////////////////////////////////////////
public:

	IntervalSlice();
	~IntervalSlice();

	// Adds an interval to this slice, combines overlapping intervals
	void addInterval(const Interval& interval);
	bool isEmpty() const;

	/////////////////////////////////////////////////////
private:
	std::vector<Interval> _intervals;

};

// Collection of interval slices, also keeps track of its position while slicing
class IntervalBread
{
	/////////////////////////////////////////////////////
public:

	IntervalBread(const Eigen::Vector2d& org, const Eigen::Vector2d& dir);
	~IntervalBread();

	// Add a slice and go to next slice
	void walk(double distance, const IntervalSlice& slice);

	// Current location
	Eigen::Vector2d getCurrent() const;

	// Starting location
	Eigen::Vector2d getAnchor() const;

	// 2D rotation matrix
	Eigen::Matrix2d getRot() const;

	// Goes through all intervals in a logical order (with as little jumps as possible), tells whether we are jumping
	void depthSearch(std::function<void(const Eigen::Vector2d&, const Eigen::Vector2d&, bool)> func, double margin, double threshold) const;

	/////////////////////////////////////////////////////
private:
	double _distance;

	Eigen::Matrix2d _rot;
	Eigen::Vector2d _org;
	std::vector<IntervalSlice> _slices;
};


// Collection of vertices that make an outline
class OutlineSlice
{
	friend class OutlineBread;

	/////////////////////////////////////////////////////
public:

	OutlineSlice(int32_t primary, int32_t secondary);
	~OutlineSlice();

	// Combines two slices if their ends match, keeps order. Returns false if ends didn't match.
	bool combine(const OutlineSlice& other);

	/////////////////////////////////////////////////////
private:
	std::deque<int32_t> _vertices;

};

// Collection of outline slices
class OutlineBread
{
	/////////////////////////////////////////////////////
public:

	OutlineBread();
	~OutlineBread();

	// Adds outline slice to this bread, tries to combine it with other slices to make a complete outline
	void addOutline(const OutlineSlice& outline);

	// Goes through all outline pieces, gives 3 vertex indices: previous, current, next
	void depthSearch(std::function<void(const Eigen::Vector3i& line, int32_t index)> func) const;

	/////////////////////////////////////////////////////
private:
	std::vector<OutlineSlice> _slices;

};



class Triangle
{
	friend class TriangleBread;

	/////////////////////////////////////////////////////
public:

	Triangle(const Eigen::Vector2d vertices[3]);
	~Triangle();

	/////////////////////////////////////////////////////
private:
	Eigen::Vector2d _vertices[3];
	Eigen::Vector2d _center;
};

class TriangleBread
{
	/////////////////////////////////////////////////////
public:

	TriangleBread();
	~TriangleBread();

	// Adds triangle to this bread
	void addTriangle(const Eigen::Vector2d vertices[3]);

	// Sorts triangles trying to minimize number of long distances between triangles
	void sortTriangles();

	// Iterates through all triangles and draws lines that go wround these triangles, generates a given amount of layers.
	void depthSearch(std::function<void(const Eigen::Vector2d&, const Eigen::Vector2d&, bool)> func, int32_t layers, double margin) const;

	/////////////////////////////////////////////////////
private:
	std::vector<Triangle> _triangles;

};