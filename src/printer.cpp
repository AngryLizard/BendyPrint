#include "printer.h"


float Printer::extruderDistance(const Vec& pos) const
{
	const float deltaX = pos.x() - extruder.pos.x();
	const float deltaY = pos.y() - extruder.pos.y();
	return sqrtf(deltaX * deltaX + deltaY * deltaY);
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

void Printer::writeMove(const Vec & pos, int32_t feed)
{
	const float dist = extruderDistance(pos);
	if (dist > EPS)
	{
		extruder.extrude += settings.lineDensity * dist;

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

void Printer::writeLayer()
{
	// Init layer
	extruder.file << "; LAYER : " << extruder.layer << std::endl;
	extruder.layer++;

	extruder.file << "M107" << std::endl;
	extruder.file << "G0 F9000 X" << extruder.pos.x() << " Y" << extruder.pos.y() << " Z" << extruder.layer * settings.layerHeight << std::endl;
	extruder.feed = 9000;
}

void Printer::writeRect(const Vec & extend, const Vec & pos)
{
	const Vec min = pos - extend;
	const Vec max = pos + extend;

	writeFast(min, true);
	writeMove(Vec(min.x(), max.y()), 600);
	writeMove(Vec(max.x(), max.y()), 600);
	writeMove(Vec(max.x(), min.y()), 600);
	writeMove(Vec(min.x(), min.y()), 600);
}


void Printer::generate(const char* filename)
{
	extruder.layer = 0;
	extruder.extrude = 0.0f;
	extruder.pos = Vec::Constant(100.0f);

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


		for (int cpy = 0; cpy < 5; cpy++)
		{
			const Vec ctr = Vec(100.0f, 50.0f + 25.0f * cpy);
			const Vec stride = settings.stride * ((float)(cpy + 1) / 5);

			// Write layers
			extruder.file << "; Layer count : " << settings.layers << std::endl;
			for (int32_t layer = 0; layer < settings.layers; layer++)
			{
				writeLayer();

				extruder.file << "; Outer " << std::endl;
				const Vec off = Vec::Constant(settings.lineWidth);
				writeRect(settings.extend - off * 0, ctr);
				writeRect(settings.extend - off * 1, ctr);

				extruder.file << "; Inner " << std::endl;
				const Vec min = ctr - (settings.extend - off * 2);
				const Vec max = ctr + (settings.extend - off * 2);
				const Vec offs = Vec::Ones() / settings.layers * settings.repetitions;
				const Vec bgn = stride.cwiseProduct(Vec(fmod(offs.x() * layer, 1.0f), fmod(offs.y() * layer, 1.0f)));
				Eigen::Vector2i totals = (max - min - bgn).cwiseQuotient(stride).cast<int32_t>();

				writeFast(min, true);
				for (int32_t i = 0; i <= totals.x(); i++)
				{
					const float x = min.x() + i * stride.x() + bgn.x();
					if (i & 1)
					{
						writeMove(Vec(x, max.y()), 780.0f);
						writeMove(Vec(x, min.y()), 780.0f);
					}
					else
					{
						writeMove(Vec(x, min.y()), 780.0f);
						writeMove(Vec(x, max.y()), 780.0f);
					}
				}

				writeFast(min, true);
				for (int32_t i = 0; i <= totals.y(); i++)
				{
					const float y = min.y() + i * stride.y() + bgn.y();
					if (i & 1)
					{
						writeMove(Vec(max.x(), y), 780.0f);
						writeMove(Vec(min.x(), y), 780.0f);
					}
					else
					{
						writeMove(Vec(min.x(), y), 780.0f);
						writeMove(Vec(max.x(), y), 780.0f);
					}
				}
			}
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
}

