#ifndef CONTAINER_h
#define CONTAINER_h

#include <vtkSmartPointer.h>
#include <memory>
#include <array>
#include <functional>
#include <optional>
#include <OpenGLRender/Model.h>
#include "Utils/Utils.h"

using Pt = std::array<double, 3>;
using Inside = std::function<bool(const Pt&)>;
using DistFunc = std::function<double(const Pt&)>;
using RadFunc = std::function<double(const double& , const double&, const double&)>;


class IContainer {
public:
	virtual ~IContainer() = default;
	virtual void render() = 0;
	virtual void gui_setup() = 0;
	virtual void create() = 0;
	virtual Bounds compute_bounds(double xDim, double yDim, double zDim) const = 0;
	virtual std::function<bool(const std::array<double, 3>&)> is_inside(
		double xDim, double yDim, double zDim) const = 0;
	virtual std::unique_ptr<DistanceEstimator> make_distance_estimator() const { return nullptr; }
	virtual std::unique_ptr<RadiusFunction> make_radius_estimator() const { return nullptr; }
};

struct BoxContainer final : IContainer {
	std::unique_ptr<BBox> model;
	float xDim{ 0.0 }, yDim{ 0.0f }, zDim{ 0.0 };
	float xMin{ 0.0f }, xMax{ 0.0f }, yMin{ 0.0f }, yMax{ 0.0f }, zMin{ 0.0f }, zMax{ 0.0f };

	void gui_setup() override {
		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along X");
		ImGui::InputFloat("x Dimensions", &xDim, 0.1f, 100.0f);

		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along Y.");
		ImGui::InputFloat("y Dimensions", &yDim, 0.1f, 100.0f);

		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along Z.");
		ImGui::InputFloat("z Dimensions", &zDim, 0.1f, 100.0f);

		xMin = 0.0;
		xMax = xDim;
		yMin = 0.0;
		yMax = yDim;
		zMin = 0.0;
		zMax = zDim;
	};

	void create() override {
		model = std::make_unique<BBox>(xMin, xMax, yMin, yMax, zMin, zMax);
	};

	void render() override {
		model->draw();
	};

	Bounds compute_bounds(double xDim, double yDim, double zDim) const override {

		double x = xDim;
		double y = yDim;
		double z = zDim;

		return{
			0.0, x, 0.0, y, 0.0, z,
			{x * 0.5f, y * 0.5f, z * 0.5f }
		};
	};

	std::function<bool(const std::array<double, 3>& pt)> is_inside(
		double xDim, double yDim, double zDim) const override {

		Bounds bounds = compute_bounds(xDim, yDim, zDim); 

		return [bounds](const std::array<double, 3>& pt) {
			return (pt[0] >= bounds.xMin && pt[0] <= bounds.xMax) &&
				(pt[1] >= bounds.yMin && pt[1] <= bounds.yMax) &&
				(pt[2] >= bounds.zMin && pt[2] <= bounds.zMax);
		};
	}

};

struct CylinderContainer final : IContainer {
	std::unique_ptr<Model> model;
	int cylinderDir{ 0 };
	float cylinderPt[3] = { 0.0f, 0.0f, 0.0f };
	float cylinderAxis[3] = { 0.0f, 0.0f, 1.0f };
	double cylinderRadius{ 2.0f };
	double cylinderHeight{ 1.0f };
	float xDim{ 0.0 }, yDim{ 0.0f }, zDim{ 0.0 };
	float xMin{ 0.0f }, xMax{ 0.0f }, yMin{ 0.0f }, yMax{ 0.0f }, zMin{ 0.0f }, zMax{ 0.0f };

	void gui_setup() override {
		ImGui::RadioButton("X", &cylinderDir, 0);
		ImGui::SameLine(); ImGui::RadioButton("Y", &cylinderDir, 1);
		ImGui::SameLine(); ImGui::RadioButton("Z", &cylinderDir, 2);
		
		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("Cylinder Radius", &cylinderRadius, 0.1, 0.0, "%.3f");

		ImGui::SetNextItemWidth(200);
		ImGui::InputDouble("Cylinder Height", &cylinderHeight, 0.1, 1.0, "%.3f");

	};

	void create() override {
		
		vtkSmartPointer<vtkCylinder> cylinder = vtkSmartPointer<vtkCylinder>::New();
		cylinder->SetRadius(cylinderRadius);

		auto cylinderAxis = unit_axis_from_dir(cylinderDir);

		cylinder->SetAxis(cylinderAxis[0], cylinderAxis[1], cylinderAxis[2]);

		double center[3]{
			 cylinderPt[0] + 0.5 * cylinderHeight * cylinderAxis[0],
			 cylinderPt[1] + 0.5 * cylinderHeight * cylinderAxis[1],
			 cylinderPt[2] + 0.5 * cylinderHeight * cylinderAxis[2],
		};

		std::cout << "Center: " << center[0] << " " << center[1] << " " << center[2] << std::endl;

		cylinder->SetCenter(center);

		model = std::make_unique<Model>(cylinder);
	};

	void render() override {
		model->draw();
	};

	Bounds compute_bounds(double xDim, double yDim, double zDim) const override {
		switch (cylinderDir) {
		case 0: // along X
			return {
				0.0, xDim, -cylinderRadius, cylinderRadius, -cylinderRadius, cylinderRadius, { xDim * 0.5, 0.0, 0.0 } };
		case 1: // along Y
			return {
				-cylinderRadius, cylinderRadius, 0.0, yDim, -cylinderRadius, cylinderRadius, { 0.0, yDim * 0.5, 0.0 } };
		default: // along Z
			return {
				-cylinderRadius, cylinderRadius, -cylinderRadius, cylinderRadius, 0.0, zDim, { 0.0, 0.0, zDim * 0.5 } };
		}
	};

	std::function<bool(const std::array<double, 3>& pt)> is_inside(
		double xDim, double yDim, double zDim) const override {

		Bounds bounds = compute_bounds(xDim, yDim, zDim);
		auto cylinderAxis = unit_axis_from_dir(cylinderDir);

		double h = cylinderHeight;
		if (h <= 0.0) {
			if (cylinderDir == 0) h = xDim;
			else if (cylinderDir == 1) h = yDim;
			else h = zDim;
		}

		const Eigen::Vector3d base(cylinderPt[0], cylinderPt[1], cylinderPt[2]);

		return [bounds, cylinderAxis, base, h, radius = cylinderRadius](const std::array<double, 3>& pt) -> bool {
			if (!is_inside_box(pt, bounds)) return false;
			return is_inside_cylinder(pt, base, cylinderAxis, radius, h);
		};
	}

};

struct MeshContainer final : IContainer {
	std::unique_ptr<Model> model;
	std::string containerFile{ "" };
	int neighbors{ 1 };
	float wallResolution{ 1.0f };
	vtkSmartPointer<vtkPolyData> containerMesh;
	double containerCenter[3]{ 0.0, 0.0, 0.0 };
	float xMin{ 0.0f }, xMax{ 0.0f }, yMin{ 0.0f }, yMax{ 0.0f }, zMin{ 0.0f }, zMax{ 0.0f };

	void gui_setup() override {

		if (ImGuiFileDialog::Instance()->Display("Select STL Container")) {
			if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
				std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
				containerFile = filePathName;
				std::cout << containerFile << std::endl;

				containerMesh = vtkSmartPointer<vtkPolyData>::New();
				vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
				reader->SetFileName(containerFile.c_str());
				reader->Update();
				containerMesh = reader->GetOutput();

				double bounds[6];
				containerMesh->GetBounds(bounds);
				xMin = bounds[0] - 1.0;
				xMax = bounds[1] + 1.0;
				yMin = bounds[2] - 1.0;
				yMax = bounds[3] + 1.0;
				zMin = bounds[4] - 1.0;
				zMax = bounds[5] + 1.0;

				// set also the mesh center as the camera pos
				float xc = (xMax + xMin) * 0.5f;
				float yc = (yMax + yMin) * 0.5f;
				float zc = (zMax + zMin) * 0.5f;
				containerCenter[0] = xc;
				containerCenter[1] = yc;
				containerCenter[2] = zc;

			}
			ImGuiFileDialog::Instance()->Close();
		}
		select_file_button("Select Geometry File", "../data/", "Select STL Container", ".stl");
		ImGui::SetNextItemWidth(200);
		ImGui::InputInt("Neighbors", &neighbors, 1, 100);
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Wall Resolution", &wallResolution, 0.1f, 5.0f);
	};

	void render() override {
		model->draw();
	};

	Bounds compute_bounds(double xDim, double yDim, double zDim) const override {

		double b[6] = { 0,0,0,0,0,0 };
		if (containerMesh) containerMesh->GetBounds(b);
		const double cx = 0.5 * (b[0] + b[1]);
		const double cy = 0.5 * (b[2] + b[3]);
		const double cz = 0.5 * (b[4] + b[5]);
		return { b[0], b[1], b[2], b[3], b[4], b[5], { cx, cy, cz } };
	};

	std::function<bool(const std::array<double, 3>& pt)> is_inside(
		double xDim, double yDim, double zDim) const override {

		Bounds bounds = compute_bounds(xDim, yDim, zDim);
		
		auto cm = containerMesh;

		return[cm](const std::array<double, 3>& pt) -> bool {
			if (!is_inside_box(pt, bounds)) return false;
			Eigen::Vector3d p(pt[0], pt[1], pt[2]);
			return is_inside_mesh(cm, p);
		};
	};

};

//@brief create an adapter for usage with the seed generators to keep it clean
struct ContainerAdapter {
	const IContainer& c;
	double xDim, yDim, zDim;
	
	Bounds get_bounds() const { return c.compute_bounds(xDim, yDim, zDim); };

	// @brief get the corresponding is_inside function
	Inside inside() const {
		return c.is_inside(xDim, yDim, zDim);
	}

	// @brief use this inside the Poisson 3D sampling to test if a point is inside
	bool is_inside(const Pt& pt) const {
		return inside()(pt);
	}

	// @brief get the distance function
	std::optional<DistFunc> distance_from() const {
		auto dEstimator = c.make_distance_estimator();
		if (dEstimator) {
			auto sp = std::shared_ptr<DistanceEstimator>(std::move(dEstimator));
			return [sp](const std::array<double, 3>& pt) {
				return std::abs(sp->compute_distance(pt));
			};
		}
		return std::nullopt;
	}

	// @brief use this as to test the distance of a point from the reference (e.g. boundaries
	// or plane)
	double get_distance(const Pt& p) const {
		if (auto f = distance_from()) return (*f)(p);
		return 0.0;
	}

};

#endif