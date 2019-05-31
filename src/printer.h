
#pragma once

#include <Eigen/Core>
#include <iostream>
#include <fstream>

#include "meshbuilder.h"

using TravelFunc = std::function<void(int32_t layer, const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)>;
using LayerFunc = std::function<void(int32_t layer)>;

using Vec = Eigen::Vector2d;

class Printer
{
public:

	struct Settings
	{
		float layerHeight;
		float lineDensity;
		float lineWidth;
		float stride;
		float threshold;

		float depth;
		float nipplescale;
		int32_t fill;
		int32_t layers;

		int32_t copies;
		float margin;

		bool fillEnabled;
		bool outlineSliceEnabled;
	} settings;

private:

	struct Extruder
	{
		std::ofstream file;
		float layer;
		float extrude;
		int32_t feed;
		Vec pos;
	} extruder;

	float extruderDistance(const Vec& pos) const;

	void writeRetract(bool lift);
	void writeExtract(bool lift);

	void writeFast(const Vec& pos, bool lift);
	void writeMove(const Vec& pos, int32_t feed, double density);

	void writeLayer(double height);
	void writeRect(const Vec& extend, const Vec& pos);

public:
	
	void travel(const Slicer* slicer, TravelFunc moveOutline, TravelFunc moveRaft, TravelFunc moveFill, LayerFunc raft = [](int32_t) {}, LayerFunc ascend = [](int32_t) {}, LayerFunc fill = [](int32_t) {});

	void generate(const Slicer* slicer, const char* filename);
};
