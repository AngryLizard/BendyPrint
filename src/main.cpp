#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>

#include "printer.h"
#include "element.h"
#include "meshbuilder.h"

void generateQuad(const Eigen::Vector2i& dim, Eigen::MatrixXd& V, Eigen::MatrixXi& F, double bend)
{
	// Compute vertices
	V = Eigen::MatrixXd::Zero((dim + Eigen::Vector2i::Ones()).prod(), 3);
	for (int32_t y = 0; y <= dim.y(); y++)
	{
		for (int32_t x = 0; x <= dim.x(); x++)
		{
			const int32_t index = y * (dim.x() + 1) + x;
			const Eigen::Vector2d p((double)x - 0.5 * dim.x(), (double)y - 0.5 * dim.y());
			V.row(index) = Vec3(p.x(), p.squaredNorm() * bend, p.y());
		}
	}

	// Compute faces
	F = Eigen::MatrixXi::Zero(dim.prod() * 2, 3);
	for (int32_t y = 0; y < dim.y(); y++)
	{
		for (int32_t x = 0; x < dim.x(); x++)
		{
			const int32_t index = (y * (dim.x() + 1) + x);
			const int32_t tmp = ((y+1) * (dim.x() + 1) + x);

			const int32_t face = (y * dim.x() + x) * 2;
			if ((x & 1) != (y & 1))
			{
				F.row(face) = Eigen::Vector3i(index, tmp + 1, tmp);
				F.row(face + 1) = Eigen::Vector3i(index, index + 1, tmp + 1);
			}
			else
			{
				F.row(face) = Eigen::Vector3i(index, index + 1, tmp);
				F.row(face + 1) = Eigen::Vector3i(tmp, index + 1, tmp + 1);
			}
		}
	}
}


int main(int argc, char *argv[])
{
	// Init the viewer
	igl::opengl::glfw::Viewer viewer;
	viewer.core.lighting_factor = 0.3f;

	// Generate mesh
	Eigen::MatrixXd pV;
	Eigen::MatrixXi pF;
	//generateQuad(Eigen::Vector2i(4, 8), pV, pF, 0.1);
	igl::readOBJ("../asset/sphere.obj", pV, pF);

	// Build shell
	Meshbuilder builder(pV, pF, 0.5);
	builder.computeElements();
	builder.renderShell(viewer);

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Settings for grid
	Printer printer;
	printer.settings.layerHeight = 0.075f;
	printer.settings.lineDensity = 0.20f;
	printer.settings.lineWidth = 0.4f;

	printer.settings.layers = 2;
	printer.settings.extend = Vec(60.0f, 10.0f);
	printer.settings.stride = Vec(4.0f, 4.0f);
	printer.settings.repetitions = 1;

	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		// Add new group
		if (ImGui::CollapsingHeader("Settings", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::DragFloat("Layer Height", &printer.settings.layerHeight, 1.0f, 0.025f, 0.1f);
			ImGui::DragFloat("Line Density", &printer.settings.lineDensity, 1.0f, 0.05f, 0.5f);
			ImGui::DragFloat("Line Width", &printer.settings.lineWidth, 1.0f, 0.01f, 1.0f);

			ImGui::DragFloat("Extend X", &printer.settings.extend.x(), 1.0f, 1.0, 100.0f);
			ImGui::DragFloat("Extend y", &printer.settings.extend.y(), 1.0f, 1.0, 100.0f);

			ImGui::DragInt("Layers", &printer.settings.layers, 1.0f, 1, 5);
			ImGui::DragInt("Repetitions", &printer.settings.repetitions, 1.0f, 1, 5);
		}

		if (ImGui::Button("Test"))
		{
			builder.gradientTest(0.0001);
		}

		if (ImGui::Button("Iterate"))
		{
			builder.gradientDescent(0.1, 10'000);
			builder.renderShell(viewer);
		}

		if (ImGui::Button("Unwrap"))
		{
			builder.computePlane(false, -3);
			builder.renderPlane(viewer);
		}

		if (ImGui::Button("Generate G-Code"))
		{
			printer.generate("generated.gcode");
		}
	};

	viewer.launch();
}
