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
	virtual void set_scale(float newScale) = 0;
	//virtual float sdf(const Vec3&) const = 0;
	std::string name = "";
	std::array<float, 4> color = {0.5, 0.5f, 0.5f , 1.0f};
	bool hidden = true;
	std::shared_ptr<SDF> sdf;
protected:
	ObjectType type = ObjectType::NoneType;
};

// ===================================================================================================
// Box Container Class
// ===================================================================================================
class BoxContainer final : public IContainer {

public:

	BoxContainer() {};
	~BoxContainer() {};

	// Updated constructor using size and origin
	BoxContainer(Vec3 size, Vec3 origin, const bool renderModeTrue = true) :
		size(size), origin(origin), renderModeTrue(renderModeTrue) {

		Bounds bds = compute_bounds();
		sdf = std::make_shared<BoxSDF>(bds);

		if (renderModeTrue) {
			create();
		}
	};

	std::unique_ptr<BBox> model;

	ObjectType get_type() const override {
		return ObjectType::BoxContainerType;
	}

	void gui_setup() override {
		bool dimensionsChanged = false;

		ImGui::ColorEdit4("Color", (float*)&color);

		ImGui::SeparatorText("Dimensions");

		ImGui::SetNextItemWidth(200);
		dimensionsChanged |= ImGui::InputFloat("Width (X) mm", &size.x, 0.1f, 100.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		dimensionsChanged |= ImGui::InputFloat("Height (Y) mm", &size.y, 0.1f, 100.0f, "%.3f");

		ImGui::SetNextItemWidth(200);
		dimensionsChanged |= ImGui::InputFloat("Depth (Z) mm", &size.z, 0.1f, 100.0f, "%.3f");

		ImGui::SeparatorText("Position");

		dimensionsChanged |= ImGui::InputFloat3("Center", (float*)&origin);
		ImGui::SetItemTooltip("Center of the Box Container");

		if (dimensionsChanged || ImGui::Button("Update")) {
			Bounds bds = compute_bounds();
			sdf = std::make_shared<BoxSDF>(bds);
			create();
		}
	};

	void create() override {

		model = std::make_unique<BBox>(size, origin);
	};

	void render() override {
		if (!model) {
			create();
		}
		if (model) {
			model->draw();
		}
	};

	Bounds compute_bounds() const override {
		float halfX = size.x * 0.5f;
		float halfY = size.y * 0.5f;
		float halfZ = size.z * 0.5f;

		return {
			origin.x - halfX, origin.x + halfX, // xMin, xMax
			origin.y - halfY, origin.y + halfY, // yMin, yMax
			origin.z - halfZ, origin.z + halfZ, // zMin, zMax
			origin                              // center
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
		return size.x * size.y * size.z;
	};

	Vec3 size{ 10.0f, 10.0f, 10.0f };

	Vec3 origin{ 5.0f, 5.0f, 5.0f };

	void set_scale(float newScale) override {};
private:

	bool renderModeTrue = true;
};

// ===================================================================================================
//@brief class to create a vertical cylindrical container along y. The bottom is 
// fixed at point (0, 0)
// ===================================================================================================
class CylinderContainer final : public IContainer {
public:
	CylinderContainer() {};
	~CylinderContainer() {};

	CylinderContainer(const float r, const float h, const bool renderMode = true) : cylinderRadius(r), cylinderHeight(h), renderMode(renderMode){
	
		// create the sdf
		sdf = std::make_shared<CylinderSDF>(
			cylinderRadius, 
			cylinderHeight, 
			Vec3(0.0f, cylinderHeight * 0.5f, 0.0f)
		);

		if (renderMode) {
			create();
		}
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

	void set_scale(float newScale) override {};

	float cylinderRadius{ 2.0f };
	float cylinderHeight{ 10.0f };

private:
	std::unique_ptr<Cylinder> model;
	
	bool renderMode = true;

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
	AbstractContainer(const std::string& fileName, bool renderMode = true);

	//@brief constructor to rebuild a mesh container from embedded geometry
	// (deduplicated vertices + faces), used when loading a .sbproj project where
	// the mesh is stored inside the project file rather than referenced by path.
	AbstractContainer(std::vector<openstl::Vec3> verts,
		std::vector<openstl::Face> faces,
		float scale = 1.0f,
		bool renderMode = true);

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
		auto meshSdf = std::dynamic_pointer_cast<MeshSDF>(sdf);

		if (meshSdf) {
			return meshSdf->is_inside_raycast(pt);
		}
		return false; 
	};

	std::shared_ptr<const SDF> get_distance_estimator() const override {
		return sdf;
	};

	bool updated = false;

	uint32_t version = 1;

	std::string fileName = "";

	// rescale the mesh about its bounding-box center (e.g. to convert a
	// CAD export's native units into mm) and rebuild the SDF/BVH
	void set_scale(float newScale) override {
		scaleFactor = newScale;
		apply_scale();
	};

	float get_scale() const {
		return scaleFactor;
	}

	// Deduplicated geometry accessors, used to embed the mesh into a .sbproj
	const std::vector<openstl::Vec3>& get_mesh_verts() const { return meshVerts; }
	const std::vector<openstl::Face>&  get_mesh_faces() const { return meshFaces; }

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
	bool renderMode = true;

	// functions
	void generate();

	// Shared construction tail: builds GL buffers, pseudonormals and the mesh
	// SDF from meshVerts/meshFaces. Called by both the file-path and the
	// embedded-geometry constructors.
	void init_from_geometry(bool renderMode);

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


// ===================================================================
// Roi Class definition
// ===================================================================
class ROI {
public:
	ROI() {};
	~ROI() {};

	// Pass size instead of 6 loose bounds
	ROI(Vec3 size, Vec3 origin, bool render = false);

	void render_model();
	std::array<float, 6> get_bounds();
	Vec3 get_center();
	Vec3 get_size() const { return size; }
	void render_properties();

	ObjectType get_type() const {
		return ObjectType::Roi;
	}

	std::string name = "";
	std::string relatedMeshName = "";
	bool hidden = false;
	ObjectType type = ObjectType::Roi;
	std::array<float, 4> color = { 1.0f, 0.0f, 0.0f, 1.0f };

private:
	Vec3 origin{ 0.0f, 0.0f, 0.0f };
	Vec3 size{ 10.0f, 10.0f, 10.0f }; 

	std::unique_ptr<BBox> model;
	bool renderMode = false;

	void create();
};

#endif