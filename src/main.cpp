#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include <igl/unproject_onto_mesh.h>

#include "printer.h"
#include "element.h"
#include "meshbuilder.h"


enum class Mode : int32_t
{
	WRAP,
	SLICE,
	FEM,
	SELECT
};

enum class SelectionMode : int32_t
{
	FIXED,
	UNFIXED,
	FORCE,
	UNFORCE
};

float scale;
float altitude;

bool primaryButton;
bool secondaryButton;

bool sliceOutlines;

float elementShear;
float elementBulk;

Mode mode;
SelectionMode selectionMode;

igl::opengl::glfw::Viewer viewer;

Printer printer;
Eigen::MatrixXd pV;
Eigen::MatrixXi pF;

Meshbuilder* builder;

// Make sure mesh is porperly displayed
void updateMode()
{
	if (mode == Mode::FEM)
	{
		builder->renderFEM(viewer, false);
	}
	else if (mode == Mode::SELECT)
	{
		builder->renderFEM(viewer, true);
	}
	else // WRAP or SLICE
	{
		builder->unwrap(0.0);
		if (mode == Mode::WRAP)
		{
			builder->getUnwrapper()->renderPlane(viewer, printer.settings.fillEnabled);
		}

		if (mode == Mode::SLICE)
		{
			const Slicer* slicer = builder->getUnwrapper()->getSlicer();

			viewer.data().clear();
			const Eigen::Vector2d nrm = Eigen::Vector2d(0.0, 1.0);

			// Translate to 3D and add as edges
			std::vector<Eigen::RowVector3d> F;
			std::vector<Eigen::RowVector3d> T;
			std::vector<Eigen::RowVector3d> C;
			std::vector<Eigen::RowVector3d> P;
			auto move = [&](int32_t layer, const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump, const Eigen::RowVector3d& col)
			{
				const double height = altitude + layer * printer.settings.layerHeight;
				if (jump) P.push_back(slicer->project(from, height));
				F.push_back(slicer->project(from, height));
				T.push_back(slicer->project(to, height));
				C.push_back(col);
			};

			// Green outlines
			TravelFunc moveOutlines = [&](int32_t layer, const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)
			{
				move(layer, from, to, jump, Eigen::RowVector3d(0.0, 1.0, 0.0));
			};

			// Alternate colour for each jump
			int32_t raftJumps = 0;
			const Eigen::RowVector3d startCol(1.0, 0.0, 0.0);
			const Eigen::RowVector3d endCol(0.0, 0.0, 1.0);
			TravelFunc moveRaft = [&](int32_t layer, const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)
			{
				move(layer, from, to, jump, startCol + (static_cast<double>(raftJumps % 5) / 5) * (endCol - startCol));
				if (jump) raftJumps++; 
			};

			// White triangles
			TravelFunc moveTriangles = [&](int32_t layer, const Eigen::Vector2d& from, const Eigen::Vector2d& to, bool jump)
			{
				move(layer + printer.settings.layers, from, to, jump, Eigen::RowVector3d(1.0, 1.0, 1.0));
			};

			// Fill line vectors
			printer.travel(slicer, moveOutlines, moveRaft, moveTriangles);

			// Display lines
			Eigen::MatrixXd FV(F.size(), 3);
			Eigen::MatrixXd TV(T.size(), 3);
			Eigen::MatrixXd CV(C.size(), 3);
			for (int i = 0; i < F.size(); i++)
			{
				FV.row(i) = F[i];
				TV.row(i) = T[i];
				CV.row(i) = C[i];
			}
			viewer.data().add_edges(FV, TV, CV);

			Eigen::MatrixXd PV(P.size(), 3);
			for (int i = 0; i < P.size(); i++)
			{
				PV.row(i) = P[i];
			}
			viewer.data().add_points(PV, Eigen::RowVector3d(1.0, 1.0, 1.0));
		}
	}
}

// Rebuilds whole mesh
void rebuild()
{
	builder->buildFromPlane(printer.settings.depth, elementShear, elementBulk);
	updateMode();
}

// Tries to select mesh at current mouse position
bool select()
{
	if (mode == Mode::SELECT && (primaryButton || secondaryButton))
	{
		const Eigen::MatrixXd& V = viewer.data().V;
		const Eigen::MatrixXi& F = viewer.data().F;

		int fid;
		Eigen::Vector3f bc;
		// Cast a ray in the view direction starting from the mouse position
		double x = viewer.current_mouse_x;
		double y = viewer.core.viewport(3) - viewer.current_mouse_y;
		if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view,
			viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
		{
			if (selectionMode == SelectionMode::FIXED)
			{
				if (primaryButton) builder->addFixedFace(fid);
				if (secondaryButton) builder->removeFixedFace(fid);
			}
			else if (selectionMode == SelectionMode::UNFIXED)
			{
				if (primaryButton) builder->removeFixedFace(fid);
				if (secondaryButton) builder->addFixedFace(fid);
			}
			else if (selectionMode == SelectionMode::FORCE)
			{
				if (primaryButton) builder->setExternalForce(fid, Eigen::Vector3d(0.0, -5.0, 0.0));
				if (secondaryButton) builder->setExternalForce(fid, Eigen::Vector3d(0.0, 0.0, 0.0));
			}
			else if (selectionMode == SelectionMode::UNFORCE)
			{
				if (primaryButton) builder->setExternalForce(fid, Eigen::Vector3d(0.0, 0.0, 0.0));
				if (secondaryButton) builder->setExternalForce(fid, Eigen::Vector3d(0.0, -5.0, 0.0));
			}
			updateMode();
			return true;
		}
	}
	return false;
}


int main(int argc, char *argv[])
{
	/*
	Eigen::Matrix3d A;
	A << 1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0;

	auto strain = [](const Eigen::Matrix3d& X) { return ( X.transpose() * X - Eigen::Matrix3d::Identity()) / 2;  };


	const double h = 0.00001;
	Eigen::Matrix3d E = strain(A);

	Tensor fE(3, 3);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Eigen::Matrix3d hA = A;
			hA(i, j) += h;
			fE(i, j) = (strain(hA) - E) / h;
		}
	}

	std::cout << A << std::endl;
	std::cout << fE * Eigen::Matrix3d::Identity() << std::endl;
	*/


	// Init the viewer
	viewer.core.lighting_factor = 0.3f;
	viewer.data().point_size = 5.0;

	scale = 20.0;
	altitude = -40.0;
	mode = Mode::FEM;
	selectionMode = SelectionMode::FIXED;

	primaryButton = false;
	secondaryButton = false;

	elementShear = 50.0;
	elementBulk = 50.0;

	// Settings for grid
	printer.settings.depth = 1.5;

	printer.settings.layerHeight = 0.3f;
	printer.settings.lineDensity = 0.20f;
	printer.settings.lineWidth = 1.1f;

	printer.settings.stride = 1.0f;
	printer.settings.threshold = 0.5f;
	printer.settings.layers = 2;

	printer.settings.copies = 1;
	printer.settings.margin = 40.0f;

	printer.settings.nipplescale = 0.85;
	printer.settings.fill = 5;
	printer.settings.fillEnabled = false;
	printer.settings.outlineSliceEnabled = true;




	// Generate mesh
	pV = Eigen::MatrixXd(4, 3);
	pF = Eigen::MatrixXi(1, 3);

	pV <<
		0.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		1.0, 0.0, 0.0,
		0.0, 0.0, 1.0;
	pF << 0, 1, 2;
	Meshbuilder* test = new Meshbuilder(pV, pF);
	test->buildFromVolume(50, 50);
	test->setVertexLocation(2, Eigen::Vector3d(2.0, 0.0, 0.0));
	test->testSimulation();
	delete test;

	/*
	*/

	/*
	pV = Eigen::MatrixXd(5, 3);
	pF = Eigen::MatrixXi(3, 3);
	pV <<
		0.9, 0.0, -4.0,
		0.0, 0.0, 1.0,
		1.0, 0.0, 0.0,
		1.1, 0.0, -4.0,
		2.0, 0.0, 1.0;
	pF <<
		0, 1, 2,
		2, 1, 4,
		2, 4, 3;
	pV *= 20.0;
	builder = new Meshbuilder(pV, pF);
	*/

	/*
	pV = Eigen::MatrixXd(3, 3);
	pF = Eigen::MatrixXi(1, 3);
	pV <<
		0.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		1.0, 0.0, 0.0;
	pF <<
		0, 1, 2;
	pV *= 20.0;
	builder = Meshbuilder(pV, pF);
	builder->buildFromPlane(printer.settings.depth);
	*/

	/*
	Eigen::Vector2i dim(4,8);
	generateQuad(dim, pV, pF, 0.1);
	Meshbuilder builder(pV, pF);
	builder->buildFromPlane(printer.settings.depth);
	*/

	//igl::readOBJ("../asset/iso.obj", pV, pF);
	//igl::readOBJ("../asset/sphere.obj", pV, pF);
	igl::readOBJ("../asset/sphere.obj", pV, pF);
	builder = new Meshbuilder(pV * scale, pF);
	/*
	*/

	rebuild();

	// Uncomment for quick FEM testing
	//builder->addFixedFace(0);
	//builder->setExternalForce(2, Eigen::Vector3d(0.0, 20, 0.0));

	// Build shell

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	std::vector<int32_t> selection;
	
	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		if (ImGui::CollapsingHeader("View Settings", ImGuiTreeNodeFlags_DefaultOpen))
		{
			std::vector<std::string> modes = {"Wrap", "Slice", "Simulation", "Selection"};
			if (ImGui::Combo("Mode", (int*)& mode, modes))
			{
				updateMode();
			}
		}

		if (ImGui::CollapsingHeader("Mesh Settings", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::DragFloat("Import scale", &scale, 0.1f, 0.1f, 100.0f);
			if (ImGui::Button("Open mesh"))
			{
				std::string fname = igl::file_dialog_open();
				if (fname.length() != 0)
				{
					if (igl::readOBJ(fname, pV, pF))
					{
						delete builder;
						builder = new Meshbuilder(pV * scale, pF);
						rebuild();
					}
				}
			}

			ImGui::DragFloat("Shell thickness", &printer.settings.depth, 0.1f, 0.1f, 50.0f);

			if (mode == Mode::FEM)
			{
				ImGui::DragFloat("Shell shear", &elementShear, 0.1f, 0.1f, 100.0f);
				ImGui::DragFloat("Shell bulk", &elementBulk, 0.1f, 0.1f, 100.0f);
			}

			if (ImGui::Button("Build mesh")) rebuild();
		}


		if (mode == Mode::WRAP || mode == Mode::SLICE)
		{
			if (ImGui::CollapsingHeader("Unwrap Settings", ImGuiTreeNodeFlags_DefaultOpen))
			{
				if (ImGui::Checkbox("Fill Enabled", &printer.settings.fillEnabled)) updateMode();
				if (mode == Mode::SLICE)
				{
					if (printer.settings.fillEnabled)
					{
						if (ImGui::CollapsingHeader("Fill settings", ImGuiTreeNodeFlags_DefaultOpen))
						{
							ImGui::DragInt("Fill Layers", &printer.settings.fill, 0.5f, 1, 5);
							ImGui::DragFloat("Fill Scale", &printer.settings.nipplescale, 0.05f, 0.1f, 1.0f);
						}
					}
				}
			}
		}

		if (mode == Mode::SLICE)
		{
			if (ImGui::CollapsingHeader("Printing Settings", ImGuiTreeNodeFlags_DefaultOpen))
			{
				if (ImGui::CollapsingHeader("General Settings", ImGuiTreeNodeFlags_DefaultOpen))
				{
					ImGui::DragFloat("Layer Height", &printer.settings.layerHeight, 1.0f, 0.025f, 0.5f);
					ImGui::DragFloat("Line Density", &printer.settings.lineDensity, 1.0f, 0.05f, 0.5f);
					ImGui::DragFloat("Line Width", &printer.settings.lineWidth, 0.05f, 0.01f, 2.0f);
					
					ImGui::DragFloat("Copy strade", &printer.settings.margin, 0.05f, 10.0f, 200.0f);
					ImGui::DragInt("Copy count", &printer.settings.copies, 1.0f, 1, 6);
				}

				if (ImGui::CollapsingHeader("Raft settings", ImGuiTreeNodeFlags_DefaultOpen))
				{
					if (ImGui::Checkbox("Slice via outlines", &printer.settings.outlineSliceEnabled)) updateMode;

					ImGui::DragFloat("Line Stride", &printer.settings.stride, 0.05f, 0.01f, 10.0f);
					ImGui::DragFloat("Line union Threshold", &printer.settings.threshold, 0.05f, 0.1f, 1.0f);

					ImGui::DragInt("Bottom Layers", &printer.settings.layers, 1.0f, 1, 5);
				}

			}
		}

		if (ImGui::CollapsingHeader("Procedures", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::Button("Refresh"))
			{
				updateMode();
			}

			if (mode == Mode::SLICE)
			{
				if (ImGui::Button("Generate G-Code"))
				{
					const Slicer* slicer = builder->getUnwrapper()->getSlicer();
					printer.generate(slicer, "generated.gcode");
				}
			}
			else if(mode == Mode::FEM)
			{

				if (ImGui::Button("Test FEM"))
				{
					builder->testSimulation();
				}

				if (ImGui::Button("Pick vertex"))
				{
					builder->setVertexLocation(0, builder->getVertexLocation(0) + Eigen::Vector3d(0.0, -5.0, 0.0));
					updateMode();
				}

				if (ImGui::Button("Iterate"))
				{
					builder->simulate(100);
					updateMode();
				}
			}
			else if (mode == Mode::SELECT)
			{
				std::vector<std::string> modes = { "Fix vertices", "Free vertices", "Add force", "Remove force" };
				ImGui::Combo("Selection mode", (int*)& selectionMode, modes);

			}
		}
	};

	viewer.callback_mouse_down =
		[&](igl::opengl::glfw::Viewer& viewer, int button, int modifier)->bool
	{
		if (button == 0)
		{
			primaryButton = true;
		}
		if (button == 2)
		{
			secondaryButton = true;
		}
		select();
		return false;
	};

	viewer.callback_mouse_up =
		[&](igl::opengl::glfw::Viewer& viewer, int button, int modifier)->bool
	{
		if (button == 0)
		{
			primaryButton = false;
		}
		if (button == 2)
		{
			secondaryButton = false;
		}
		return false;
	};

	viewer.callback_mouse_move =
		[&](igl::opengl::glfw::Viewer& viewer, int x, int y)->bool
	{

		return select();
	};

	viewer.launch();
	delete builder;
}
