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
#include "Math/Vec.h"

using Pt = std::array<double, 3>;
using Inside = std::function<bool(const Pt&)>;

class IContainer {
public:
	virtual ~IContainer() = default;
	virtual void render() = 0;
	virtual void gui_setup() = 0;
	virtual void create() = 0;
	virtual Bounds compute_bounds() const = 0;
	virtual ObjectType get_type() const = 0;
	//virtual std::function<bool(const Vec3&)> is_inside() const = 0;
	virtual bool is_inside(const Vec3& pt) const = 0;
	virtual std::shared_ptr<const SDF> get_distance_estimator() const = 0;
	//virtual float sdf(const Vec3&) const = 0;
	std::string name = "";
	bool hidden = true;
	std::shared_ptr<SDF> sdf;
protected:
	ObjectType type = ObjectType::NoneType;
};

class BoxContainer final : public IContainer {
	
public:

	BoxContainer() {};
	~BoxContainer() {};

	// add a constructor
	BoxContainer(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax) : 
		xMin(xMin), xMax(xMax), yMin(yMin), yMax(yMax), zMin(zMin), zMax(zMax) {
		
		Bounds bds = compute_bounds();

		sdf = std::make_shared<BoxSDF>(bds);
		
		create();
	};

	std::unique_ptr<BBox> model;

	ObjectType get_type() const override {
		return ObjectType::BoxContainerType;
	}

	void gui_setup() override {
		
		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along X (mm)");
		ImGui::InputFloat("x Min (mm)", &xMin, 0.1f, 100.0f, "%.3f");
		ImGui::InputFloat("x Max (mm)", &xMax, 0.1f, 100.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along Y (mm)");
		ImGui::InputFloat("y Min (mm)", &yMin, 0.1f, 100.0f, "%.3f");
		ImGui::InputFloat("y Max (mm)", &yMax, 0.1f, 100.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		ImGui::SetItemTooltip("Dimension of Scaffold Along Z (mm)");
		ImGui::InputFloat("z Min (mm)", &zMin, 0.1f, 100.0f, "%.3f");
		ImGui::InputFloat("z Max (mm)", &zMax, 0.1f, 100.0f, "%.3f");

		if (ImGui::Button("Update")) {
			Bounds bds = compute_bounds();
			sdf = std::make_shared<BoxSDF>(bds);
			create();
		}
	};

	void create() override {
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

		return{
			xMin, xMax, yMin, yMax, zMin, zMax,
			{(xMax - xMin) * 0.5f, (yMax - yMin) * 0.5f, (zMax - zMin) * 0.5f }
		};
	};

	bool is_inside(const Vec3& pt) const override {
		if (sdf->compute_distance(pt) > 1e-4) {
			return false;
		}
		else {
			return true;
		}
	};

	virtual std::shared_ptr<const SDF> get_distance_estimator() const {
		return sdf;
	};

private:
	float xMin{ 0.0f }, xMax{ 10.0f }, yMin{ 0.0f }, yMax{ 10.0f }, zMin{ 0.0f }, zMax{ 10.0f };
};

//@brief class to create a vertical cylindrical container along y. The bottom is 
// fixed at point (0, 0)
class CylinderContainer final : public IContainer {
public:
	CylinderContainer() {};
	~CylinderContainer() {};

	CylinderContainer(const float r, const float h) :cylinderRadius(r), cylinderHeight(h){
	
		// create the sdf
		sdf = std::make_shared<CylinderSDF>(cylinderRadius, cylinderHeight, Vec3(0.0f, cylinderHeight * 0.5f, 0.0f));

		create();
	};

	void gui_setup() override {

		static bool flag = false;

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Cylinder Radius", &cylinderRadius, 0.1f, 0.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("Cylinder Height", &cylinderHeight, 0.1f, 1.0f, "%.3f");

		if (ImGui::Button("Update")) {
			
			sdf = std::make_shared<CylinderSDF>(cylinderRadius, cylinderHeight,
				Vec3(0.0f, cylinderHeight * 0.5f, 0.0f));

			create();
		}
	};

	ObjectType get_type() const override {
		return ObjectType::CylinderContainerType;
	}

	virtual std::shared_ptr<const SDF> get_distance_estimator() const override {
		return sdf;
	};

private:
	std::unique_ptr<Cylinder> model;
	
	float cylinderRadius{ 2.0f };
	float cylinderHeight{ 10.0f };	

	void create() override {

		model = std::make_unique<Cylinder>(
			Vec3{ 0.0f, 0.0f, 0.0f },
			Vec3{ 0.0f, 0.1f, 0.0f },
			cylinderHeight, cylinderRadius);
	};

	void render() override {
		if (model) {
			model->draw();
		}
	};

	Bounds compute_bounds() const override {

		Vec3 center = { 0.0f, cylinderHeight * 0.5f, 0.0f };

		return {
			-cylinderRadius, cylinderRadius,
			0.0f, cylinderHeight,
			-cylinderRadius, cylinderRadius, center};
	};

	bool is_inside(const Vec3& pt) const override {

		if (sdf->compute_distance(pt) > 1e-4) {
			return false;
		}
		else {
			return true;
		}
	};
};

//struct MeshContainer final : IContainer {
//	std::unique_ptr<Model> model;
//	std::string containerFile{ "" };
//	int neighbors{ 1 };
//	float wallResolution{ 1.0f };
//	vtkSmartPointer<vtkPolyData> containerMesh;
//	double containerCenter[3]{ 0.0, 0.0, 0.0 };
//	float xMin{ 0.0f }, xMax{ 0.0f }, yMin{ 0.0f }, yMax{ 0.0f }, zMin{ 0.0f }, zMax{ 0.0f };
//
//	ObjectType get_type() const override {
//		return ObjectType::AbstractContainerType;
//	}
//
//	void gui_setup() override {
//
//		select_file_button("Select Geometry File", "../data/", "Select STL Container", ".stl");
//		ImGui::SetNextItemWidth(200);
//		ImGui::InputInt("Neighbors", &neighbors, 1, 100);
//		ImGui::SetNextItemWidth(200);
//		ImGui::InputFloat("Wall Resolution", &wallResolution, 0.1f, 5.0f);
//
//		ImVec2 minSize = ImVec2(500.0f, 500.0f);
//		if (ImGuiFileDialog::Instance()->Display("Select STL Container", ImGuiWindowFlags_NoCollapse, minSize)) {
//			if (ImGuiFileDialog::Instance()->IsOk()) { // action if OK
//				std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
//				containerFile = filePathName;
//				std::cout << containerFile << std::endl;
//
//				containerMesh = vtkSmartPointer<vtkPolyData>::New();
//				vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
//				reader->SetFileName(containerFile.c_str());
//				reader->Update();
//				containerMesh = reader->GetOutput();
//
//				double bounds[6];
//				containerMesh->GetBounds(bounds);
//				xMin = bounds[0] - 1.0;
//				xMax = bounds[1] + 1.0;
//				yMin = bounds[2] - 1.0;
//				yMax = bounds[3] + 1.0;
//				zMin = bounds[4] - 1.0;
//				zMax = bounds[5] + 1.0;
//
//				// set also the mesh center as the camera pos
//				float xc = (xMax + xMin) * 0.5f;
//				float yc = (yMax + yMin) * 0.5f;
//				float zc = (zMax + zMin) * 0.5f;
//				containerCenter[0] = xc;
//				containerCenter[1] = yc;
//				containerCenter[2] = zc;
//
//			}
//			ImGuiFileDialog::Instance()->Close();
//		}
//	};
//
//	void create() override {
//		if (containerMesh) {
//			model = std::make_unique<Model>(containerMesh);
//		}
//	};
//
//	void render() override {
//		if (model) {
//			model->draw();
//		}
//	};
//
//	Bounds compute_bounds() const override {
//
//		double b[6] = { 0,0,0,0,0,0 };
//		if (containerMesh) containerMesh->GetBounds(b);
//		const double cx = 0.5 * (b[0] + b[1]);
//		const double cy = 0.5 * (b[2] + b[3]);
//		const double cz = 0.5 * (b[4] + b[5]);
//		return { b[0], b[1], b[2], b[3], b[4], b[5], { cx, cy, cz } };
//	};
//
//	std::function<bool(const std::array<double, 3>& pt)> is_inside() const override {
//
//		Bounds bounds = compute_bounds();
//
//		auto cm = containerMesh;
//
//		return[cm](const std::array<double, 3>& pt) -> bool {
//			Eigen::Vector3d p(pt[0], pt[1], pt[2]);
//			return is_inside_mesh(cm, p);
//		};
//	};
//
//	std::shared_ptr<const DistanceEstimator> get_distance_estimator() const override {
//		return std::make_shared<MeshDistEstimator>(containerMesh);
//	};
//};

//@brief create an adapter for usage with the seed generators to keep it clean
//struct ContainerAdapter {
//	const IContainer& c;
//	double xDim, yDim, zDim;
//
//	ContainerAdapter(const IContainer& c_, double x, double y, double z)
//		: c(c_), xDim(x), yDim(y), zDim(z) {
//	}
//
//	Bounds get_bounds() const { return c.compute_bounds(); };
//
//	// @brief get the corresponding is_inside function
//	Inside inside() const {
//		return c.is_inside();
//	}
//
//	// @brief use this inside the Poisson 3D sampling to test if a point is inside
//	bool is_inside(const Pt& pt) const {
//		return inside()(pt);
//	};
//
//	std::shared_ptr<const DistanceEstimator> get_estimator() const {
//
//		// if the estimator is not cached estimate it and cache it
//		if (!cached_estimator) {
//			cached_estimator = c.get_distance_estimator();
//		}
//		return cached_estimator;
//	};
//
//private:
//	mutable std::shared_ptr<const DistanceEstimator> cached_estimator;
//};




#endif