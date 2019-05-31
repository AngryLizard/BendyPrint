#include "printer.h"


float Printer::extruderDistance(const Vec& pos) const
{
	const float deltaX = pos.x() - extruder.pos.x();
	const float deltaY = pos.y() - extruder.pos.y();
	return sqrt(deltaX * deltaX + deltaY * deltaY);
}

void Printer::writeRetract(bool lift)
{
	extruder.feed = 2400;
	if (lift)
	{
		extruder.file << "G1 F" << extruder.feed << " E" << extruder.extrude - 3.0f << " Z" << (extruder.layer + 3) * settings.layerHeight << std::endl;
	}
	else
	{
		extruder.file << "G1 F" << extruder.feed << " E" << extruder.extrude - 3.0f << std::endl;
	}
}

void Printer::writeExtract(bool lift)
{
	extruder.feed = 2400;
	if (lift)
	{
		extruder.file << "G1 F" << extruder.feed << " E" << extruder.extrude << " Z" << (extruder.layer) * settings.layerHeight << std::endl;
	}
	else
	{
		extruder.file << "G1 F" << extruder.feed << " E" << extruder.extrude << std::endl;
	}
}

void Printer::writeFast(const Vec & pos, bool lift)
{
	if (extruderDistance(pos) > EPS)
	{
		writeRetract(lift);
		extruder.feed = 9000;
		extruder.file << "G0 F" << extruder.feed << " X" << pos.x() << " Y" << pos.y() << std::endl;
		extruder.pos = pos;
		writeExtract(lift);
	}
}

void Printer::writeMove(const Vec & pos, int32_t feed, double density)
{
	const float dist = extruderDistance(pos);
	if (dist > EPS)
	{
		extruder.extrude += settings.lineDensity * dist * density;

		if (feed != extruder.feed)
		{
			extruder.file << "G1 F" << feed << " X" << pos.x() << " Y" << pos.y() << " E" << extruder.extrude << std::endl;
			extruder.feed = feed;
		}
		else
		{
			extruder.file << "G1 X" << pos.x() << " Y" << pos.y() << " E" << extruder.extrude << std::endl;
		}
		extruder.pos = pos;
	}
}

void Printer::writeLayer(double height)
{
	// Init layer
	extruder.file << "; LAYER : " << extruder.layer << std::endl;
	extruder.layer += settings.layerHeight * height;

	extruder.file << "M107" << std::endl;
	extruder.file << "G0 F9000 X" << extruder.pos.x() << " Y" << extruder.pos.y() << " Z" << extruder.layer << std::endl;
	extruder.feed = 9000;
}

void Printer::writeRect(const Vec & extend, const Vec & pos)
{
	const Vec min = pos - extend;
	const Vec max = pos + extend;

	writeFast(min, true);
	writeMove(Vec(min.x(), max.y()), 600, 1.0);
	writeMove(Vec(max.x(), max.y()), 600, 1.0);
	writeMove(Vec(max.x(), min.y()), 600, 1.0);
	writeMove(Vec(min.x(), min.y()), 600, 1.0);
}

void Printer::travel(const Slicer* slicer, TravelFunc moveOutline, TravelFunc moveRaft, TravelFunc moveFill, LayerFunc raft, LayerFunc ascend, LayerFunc fill)
{
	// Travel layers
	for (int32_t layer = 0; layer < settings.layers; layer++)
	{
		ascend(layer);

		const OutlineBread outlineBread = slicer->createOutlineBread();
		slicer->outlineSearch(outlineBread, [&](const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump) {
			moveOutline(layer, from, to, jump);
			}, settings.lineWidth * 0.6);

		raft(layer);
		const double angle = M_PI / (settings.layers) * layer;

		const Eigen::Vector2d nrm = Eigen::Vector2d(std::cos(angle), std::sin(angle));
		if (settings.outlineSliceEnabled)
		{
			const IntervalBread intervalBread = slicer->createIntervalBread(nrm.normalized(), outlineBread, settings.stride);
			auto move = [&](const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump) { moveRaft(layer, from, to, jump); };
			intervalBread.depthSearch(move, 0.0, settings.threshold);
		}
		else
		{
			const IntervalBread intervalBread = slicer->createIntervalBread(nrm.normalized(), settings.stride);
			auto move = [&](const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump) { moveRaft(layer, from, to, jump); };
			intervalBread.depthSearch(move, 0.0, settings.threshold);
		}
	}

	if (settings.fillEnabled)
	{
		const int32_t layers = (int32_t)(settings.depth / settings.layerHeight);
		for (int32_t layer = 0; layer < layers; layer++)
		{
			fill(layer);

			const double height = settings.layerHeight * layer;
			const TriangleBread triangleBread = slicer->createTriangleBread(settings.nipplescale, height);
			auto move = [&](const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump) { moveFill(layer, from, to, jump); };
			triangleBread.depthSearch(move, settings.fill, settings.lineWidth);
		}
	}
}


void Printer::generate(const Slicer* slicer, const char* filename)
{
	Settings memory = settings;
	Vec ctr = Vec(100.0f, 100.0f);

	const int32_t speed = 960;

	extruder.layer = 0.0f;
	extruder.extrude = 0.0f;
	extruder.pos = ctr;

	const std::ios_base::openmode mode = std::ios::out | std::ios::binary | std::ios::trunc;
	extruder.file = std::ofstream(filename, mode);

	if (extruder.file.is_open())
	{
		// Prologue template
		std::ifstream prologue("../templates/Prologue.txt", std::ios::in);
		if (prologue.is_open())
		{
			extruder.file << prologue.rdbuf();
		}
		else
		{
			std::cerr << "Opening prologue failed!" << std::endl;
		}
		prologue.close();
		
		ctr.x() -= 0.5f * settings.margin * ((double)settings.copies - 0.5);
		for (int copy = 0; copy < settings.copies; copy++)
		{
			if (settings.copies > 1)
			{
				// Vary some setting for each copy
				const double variation = 0.5 + 1.0 * ((double)copy / (settings.copies - 1));
				settings.lineDensity = memory.lineDensity * variation;
			}

			// Extrudes or moves on a jump
			TravelFunc move = [&](int32_t layer, const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)
			{
				const Eigen::Vector2d start = ctr + from;
				const Eigen::Vector2d end = ctr + to;
				if (jump)
				{
					writeFast(start, true);
				}
				else
				{
					writeMove(start, speed, 1.0);
				}
				writeMove(end, speed, 1.0);
			};

			// Change direction unformly
			LayerFunc raft = [&](int32_t layer)
			{
				extruder.file << "; Inner " << std::endl;
			};

			// Write half a layer for the raft
			LayerFunc ascend = [&](int32_t layer)
			{
				writeLayer(0.25f);// layer == 0 ? 0.5f : 1.0f);
				extruder.file << "; Outer " << std::endl;
			};

			// Write layer
			LayerFunc fill = [&](int32_t layer)
			{
				writeLayer(2.0);
			};

			// Write layers
			extruder.file << "; Layer count : " << settings.layers << std::endl;
			travel(slicer, move, move, move, raft, ascend, fill);


			ctr.x() += settings.margin;
		}


		// Ending template
		std::ifstream ending("../templates/Ending.txt", std::ios::in);
		if (ending.is_open())
		{
			extruder.file << ending.rdbuf();
		}
		else
		{
			std::cerr << "Opening ending failed!" << std::endl;
		}
		ending.close();

		extruder.file.close();
		std::cerr << "Writing extruder file finished!" << std::endl;
	}
	else
	{
		std::cerr << "Opening file " << filename << " failed!" << std::endl;
	}
	settings = memory;
}

