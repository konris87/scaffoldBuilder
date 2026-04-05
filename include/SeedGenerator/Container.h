#ifndef CONTAINER_h
#define CONTAINER_h

#include <memory>
#include <array>
#include <string>
#include <functional>
#include <optional>
#include "Utils/Utils.h"
#include "DistanceCalculator.h"
#include "RadiusCalculator.h"
#include "Math/Vec.h"
#include <Openstl/core/stl.h>
#include <OpenGLRender/Model.h>

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
	virtual float get_volume() const = 0;
	//virtual float sdf(const Vec3&) const = 0;
	std::string name = "";
	std::array<float, 4> color = {0.5, 0.5f, 0.5f , 1.0f};
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
		ImGui::InputFloat("x Min (mm)", &xMin, 0.1f, 100.0f, "%.3f");
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("x Max (mm)", &xMax, 0.1f, 100.0f, "%.3f");
		ImGui::SetItemTooltip("Dimension of Scaffold Along X (mm)");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("y Min (mm)", &yMin, 0.1f, 100.0f, "%.3f");
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("y Max (mm)", &yMax, 0.1f, 100.0f, "%.3f");
		ImGui::SetItemTooltip("Dimension of Scaffold Along Y (mm)");

		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("z Min (mm)", &zMin, 0.1f, 100.0f, "%.3f");
		ImGui::SetNextItemWidth(200);
		ImGui::InputFloat("z Max (mm)", &zMax, 0.1f, 100.0f, "%.3f");
		ImGui::SetItemTooltip("Dimension of Scaffold Along Z (mm)");

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

	std::shared_ptr<const SDF> get_distance_estimator() const override {
		return sdf;
	};

	float get_volume() const override {
	
		return (xMax - xMin) * (yMax - yMin) * (zMax - zMin);
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

	CylinderContainer(const float r, const float h) : cylinderRadius(r), cylinderHeight(h){
	
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

	std::shared_ptr<const SDF> get_distance_estimator() const override {
		return sdf;
	};

	float get_volume() const override {
		const float PI = 3.14159265359f;

		return PI * (cylinderRadius * cylinderRadius) * cylinderHeight;
	};

	bool is_inside(const Vec3& pt) const override {

		if (sdf->compute_distance(pt) >= 1e-4) {
			return false;
		}
		else {
			return true;
		}
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
};

class AbstractContainer : public IContainer {
public:
	AbstractContainer() {};
	~AbstractContainer() {};

	//@brief constructor to load an stl mesh using opensstl
	AbstractContainer(const std::string& fileName);

	void gui_setup() override;

	ObjectType get_type() const override {
		return ObjectType::AbstractContainerType;
	}

	float get_volume() const override {
		// we can use the tetrahedron formula to get the volume
		float volume = 0.0f;
		float signedVolume = 0.0f;
	
		for (const auto& tri : triangles) {

			// grab the three vertices
			const Vec3& v1 = tri.v1;
			const Vec3& v2 = tri.v2;
			const Vec3& v3 = tri.v3;

			// estimate the surface area
			Vec3 edge1 = Vec3(v2.x, v2.y, v2.z) - Vec3(v1.x, v1.y, v1.z);
			Vec3 edge2 = Vec3(v3.x, v3.y, v3.z) - Vec3(v1.x, v1.y, v1.z);


			// get the volume using the signed tetrahedron 
			float vol = Vec3(v1.x, v1.y, v1.z).dot(Vec3(v2.x, v2.y, v2.z).cross(Vec3(v3.x, v3.y, v3.z)));
			signedVolume += vol;
		}

		// we need to divide the volume by 1/6
		volume = std::abs(signedVolume) / 6.0f;
		return volume;
	};

	void render() override {
		if (model) {
			model->draw();
		}
	};

	bool is_inside(const Vec3& pt) const override {

		if (sdf->compute_distance(pt) >= 1e-4) {
			return false;
		}
		else {
			return true;
		}
	};

	std::shared_ptr<const SDF> get_distance_estimator() const override {
		return sdf;
	};
	bool updated = false;

	uint32_t version = 1;

private:
	std::array<float, 6> bounds = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

	// these are for opengl
	std::vector<float> vertices;
	std::vector<unsigned int> indices;
	std::vector<float> vertexNormals;
	std::vector<unsigned int> edgeIndices;
	std::vector<Triangle> triangles;
	std::set<std::pair<unsigned int, unsigned int>> edgeSet;
	std::vector<std::vector<unsigned int>> adjacency;
	std::vector<std::vector<unsigned int>> faceAdjacency;
	std::vector<Vec3> pseudonormals;
	std::unique_ptr<Model> model;
	std::vector<openstl::Vec3> meshVerts;
	std::vector<openstl::Face> meshFaces;
	float scaleFactor = 1.0f;

	// functions
	void generate();

	void create() override {

		model = std::make_unique<Model>(vertices, indices, vertexNormals);
	};

	Bounds compute_bounds() const override {
		
		// find the box center 
		Vec3 center = {
			(bounds[1] + bounds[0]) * 0.5f,
			(bounds[3] + bounds[2]) * 0.5f,
			(bounds[5] + bounds[4]) * 0.5f };

		return { bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5], center };
	};
	
	void estimate_pseudonormals();

	void apply_scale();
};


#endif