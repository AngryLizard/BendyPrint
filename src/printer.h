
#pragma once

#include <Eigen/Core>
#include <iostream>
#include <fstream>

#define EPS 0.0001f

using Vec = Eigen::Vector2f;


class Printer
{
public:

	struct Settings
	{
		float layerHeight;
		float lineDensity;
		float lineWidth;

		int32_t layers;
		Vec extend;
		Vec stride;
		int32_t repetitions;
	} settings;

private:

	struct Extruder
	{
		std::ofstream file;
		int32_t layer;
		float extrude;
		int32_t feed;
		Vec pos;
	} extruder;

	float extruderDistance(const Vec& pos) const;

	void writeRetract(bool lift);
	void writeExtract(bool lift);

	void writeFast(const Vec& pos, bool lift);
	void writeMove(const Vec& pos, int32_t feed);

	void writeLayer();
	void writeRect(const Vec& extend, const Vec& pos);

public:

	void generate(const char* filename);
};
