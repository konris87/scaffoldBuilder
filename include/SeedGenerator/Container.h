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
			{(xMax + xMin) * 0.5f, (yMax + yMin) * 0.5f, (zMax + zMin) * 0.5f }
		};
	};

	bool is_inside(const Vec3& pt) const override {
		if (sdf->compute_distance(pt) >= 1e-4) {
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

		if (sdf->compute_distance(pt) >= 1e-4) {
			return false;
		}
		else {
			return true;
		}
	};
};


#endif