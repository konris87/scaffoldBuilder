#ifndef CONTAINER_h
#define CONTAINER_h

#include <memory>
#include <array>
#include <functional>
#include <optional>
#include <OpenGLRender/Model.h>
#include "Utils/Utils.h"
#include "DistanceCalculator.h"
#include "RadiusCalculator.h"
#include <vtkImplicitFunction.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkBox.h>
#include <vtkCylinder.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>

using Pt = std::array<double, 3>;
using Inside = std::function<bool(const Pt&)>;

class IContainer {
public:
	virtual ~IContainer() = default;
	virtual void render() = 0;
	virtual void gui_setup() = 0;
	virtual void create() = 0;
	virtual Bounds compute_bounds() const = 0;
	virtual std::function<bool(const std::array<double, 3>&)> is_inside() const = 0;
	virtual std::shared_ptr<const DistanceEstimator> get_distance_estimator() const = 0;
};

struct BoxContainer final : IContainer {
	std::unique_ptr<BBox> model;
	float xDim{ 10.0 }, yDim{ 10.0f }, zDim{ 10.0 };
	float xMin{ 0.0f }, xMax{ 0.0f }, yMin{ 0.0f }, yMax{ 0.0f }, zMin{ 0.0f }, zMax{ 0.0f };

	// add a constructor
	BoxContainer() = default;
	//BoxContainer(float xDim, float yDim, float zDim) : xDim(xDim), yDim(yDim), zDim(zDim) {
	//	create();
	//};

	void gui_setup() override {
		
		bool flag = false;
		
		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along X");
		flag |= ImGui::InputFloat("x Dimensions", &xDim, 0.1f, 100.0f);

		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along Y.");
		flag |= ImGui::InputFloat("y Dimensions", &yDim, 0.1f, 100.0f);

		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along Z.");
		flag |= ImGui::InputFloat("z Dimensions", &zDim, 0.1f, 100.0f);

		if (flag) {
			create();
			flag = false;
		}
	};

	void create() override {
		xMin = 0.0;
		xMax = xDim;
		yMin = 0.0;
		yMax = yDim;
		zMin = 0.0;
		zMax = zDim;
		model = std::make_unique<BBox>(xMin, xMax, yMin, yMax, zMin, zMax);
	};

	void render() override {
		if (!model){
			create();
		}
		if (model) {
			model->draw();
		}
	};

	Bounds compute_bounds() const override {

		double x = xDim;
		double y = yDim;
		double z = zDim;

		return{
			0.0, x, 0.0, y, 0.0, z,
			{x * 0.5f, y * 0.5f, z * 0.5f }
		};
	};

	std::function<bool(const std::array<double, 3>& pt)> is_inside() const override {

		Bounds bounds = compute_bounds(); 

		return [bounds](const std::array<double, 3>& pt) {
			return (pt[0] >= bounds.xMin && pt[0] <= bounds.xMax) &&
				(pt[1] >= bounds.yMin && pt[1] <= bounds.yMax) &&
				(pt[2] >= bounds.zMin && pt[2] <= bounds.zMax);
		};
	}

	std::shared_ptr<const DistanceEstimator> get_distance_estimator() const override {
	
		vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
		box->SetBounds(xMin, xMax, yMin, yMax, zMin, zMax);

		return std::make_shared<ImplicitFunctionDistEstimator>(box);
	};

};

struct CylinderContainer final : IContainer {
	std::unique_ptr<Cylinder> model;
	int cylinderDir{ 0 };
	float cylinderPt[3] = { 0.0f, 0.0f, 0.0f };
	float cylinderAxis[3] = { 0.0f, 0.0f, 1.0f };
	double cylinderRadius{ 2.0f };
	double cylinderHeight{ 1.0f };
	float xDim{ 0.0 }, yDim{ 0.0f }, zDim{ 0.0 };
	float xMin{ 0.0f }, xMax{ 0.0f }, yMin{ 0.0f }, yMax{ 0.0f }, zMin{ 0.0f }, zMax{ 0.0f };

	void gui_setup() override {

		static bool flag = false;

		flag |= ImGui::RadioButton("X", &cylinderDir, 0);
		ImGui::SameLine(); 
		flag |= ImGui::RadioButton("Y", &cylinderDir, 1);
		ImGui::SameLine();
		flag |= ImGui::RadioButton("Z", &cylinderDir, 2);
		
		ImGui::SetNextItemWidth(200);
		flag |= ImGui::InputDouble("Cylinder Radius", &cylinderRadius, 0.1, 0.0, "%.3f");

		ImGui::SetNextItemWidth(200);
		flag |= ImGui::InputDouble("Cylinder Height", &cylinderHeight, 0.1, 1.0, "%.3f");

		switch (cylinderDir) {
			case 0: {
				xDim = cylinderHeight;
				yDim = cylinderRadius;
				zDim = cylinderRadius;
				break;
				}
			case 1: {
				xDim = cylinderRadius;
				yDim = cylinderHeight;
				zDim = cylinderRadius;
				break;
				}
			case 2: {
				xDim = cylinderRadius;
				yDim = cylinderRadius;
				zDim = cylinderHeight;
				break;
			}
		}

		if (flag) {
			create();
			flag = false;
		}
	};

	void create() override {

		//vtkSmartPointer<vtkCylinder> cylinder = vtkSmartPointer<vtkCylinder>::New();
		//cylinder->SetRadius(cylinderRadius);

		auto cylinderAxis = unit_axis_from_dir(cylinderDir);

		model = std::make_unique<Cylinder>(
			Vec3{ cylinderPt[0], cylinderPt[1], cylinderPt[2] },
			Vec3{ (float)cylinderAxis[0], (float)cylinderAxis[1], (float)cylinderAxis[2] },
			cylinderHeight, cylinderRadius);
	};

	void render() override {
		if (model) {
			model->draw();
		}
	};

	Bounds compute_bounds() const override {
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

	std::function<bool(const std::array<double, 3>& pt)> is_inside() const override {

		Bounds bounds = compute_bounds();
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

	std::shared_ptr<const DistanceEstimator> get_distance_estimator() const override {

		vtkSmartPointer<vtkCylinder> cylinder = vtkSmartPointer<vtkCylinder>::New();
		cylinder->SetRadius(cylinderRadius);

		auto cylinderAxis = unit_axis_from_dir(cylinderDir);

		cylinder->SetAxis(cylinderAxis[0], cylinderAxis[1], cylinderAxis[2]);

		double center[3]{
			 cylinderPt[0] + 0.5 * cylinderHeight * cylinderAxis[0],
			 cylinderPt[1] + 0.5 * cylinderHeight * cylinderAxis[1],
			 cylinderPt[2] + 0.5 * cylinderHeight * cylinderAxis[2],
		};

		cylinder->SetCenter(center);

		return std::make_shared<ImplicitFunctionDistEstimator>(cylinder);
	};

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

		select_file_button("Select Geometry File", "../data/", "Select STL Container", ".stl");
		ImGui::SetNextItemWidth(200);
		ImGui::InputInt("Neighbors", &neighbors, 1, 100);
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Wall Resolution", &wallResolution, 0.1f, 5.0f);

		ImVec2 minSize = ImVec2(500.0f, 500.0f);
		if (ImGuiFileDialog::Instance()->Display("Select STL Container", ImGuiWindowFlags_NoCollapse, minSize)) {
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
	};

	void create() override {
		if (containerMesh) {
			model = std::make_unique<Model>(containerMesh);
		}
	};

	void render() override {
		if (model) {
			model->draw();
		}
	};

	Bounds compute_bounds() const override {

		double b[6] = { 0,0,0,0,0,0 };
		if (containerMesh) containerMesh->GetBounds(b);
		const double cx = 0.5 * (b[0] + b[1]);
		const double cy = 0.5 * (b[2] + b[3]);
		const double cz = 0.5 * (b[4] + b[5]);
		return { b[0], b[1], b[2], b[3], b[4], b[5], { cx, cy, cz } };
	};

	std::function<bool(const std::array<double, 3>& pt)> is_inside() const override {

		Bounds bounds = compute_bounds();
		
		auto cm = containerMesh;

		return[cm](const std::array<double, 3>& pt) -> bool {
			Eigen::Vector3d p(pt[0], pt[1], pt[2]);
			return is_inside_mesh(cm, p);
		};
	};

	std::shared_ptr<const DistanceEstimator> get_distance_estimator() const override {
		return std::make_shared<MeshDistEstimator>(containerMesh);
	};
};

//@brief create an adapter for usage with the seed generators to keep it clean
struct ContainerAdapter {
	const IContainer& c;
	double xDim, yDim, zDim;

	ContainerAdapter(const IContainer& c_, double x, double y, double z)
		: c(c_), xDim(x), yDim(y), zDim(z) {
	}

	Bounds get_bounds() const { return c.compute_bounds(); };

	// @brief get the corresponding is_inside function
	Inside inside() const {
		return c.is_inside();
	}

	// @brief use this inside the Poisson 3D sampling to test if a point is inside
	bool is_inside(const Pt& pt) const {
		return inside()(pt);
	};

	std::shared_ptr<const DistanceEstimator> get_estimator() const {

		// if the estimator is not cached estimate it and cache it
		if (!cached_estimator) {
			cached_estimator = c.get_distance_estimator();
		}
		return cached_estimator;
	};

private:
	mutable std::shared_ptr<const DistanceEstimator> cached_estimator;
};

#endif