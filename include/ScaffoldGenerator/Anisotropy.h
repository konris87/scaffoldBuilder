#ifndef ANISOTROPY_H
#define ANISOTROPY_H

/* 
Implement anisotropy enforcement 
======================================

*/

#include <vector>
#include "Eigen/Dense"
#include "Math/Vec.h"
#include "Utils/Utils.h"
#include "OpenGlRender/Model.h"
#include "Misc/Imgui_Stdlib.h"

// ===================================================================
// Anisotropy Source Class definition
// ===================================================================

class AnisotropySource {

public:
    AnisotropySource(
        Vec3 center = {0.0f, 0.0f, 0.0f},
        Vec3 direction = {1.0f, 0.0f, 0.0f},
        float angle = 0.0f,
        Vec3 stretch = {1.0f, 1.0f, 1.0f},
        float sigma = 1.1f
    ){};

    Vec3 origin; // set initially to zero
    Vec3 direction = {1.0f, 0.0f, 0.0f};
    float angle = 0.0f;
    Vec3 stretch = {1.0f, 1.0f, 1.0f};
    float sigma = 1.1f; // gaussian falloff parameter
    Eigen::Matrix3f M;

    std::string name = "";
	float color[4] = {0.0f, 1.0f, 0.0f, 0.4f};
	float lineColor[4] = {1.0f, 0.0f, 0.0f, 0.9f};
    bool hidden = false;

    void render_sphere_model() {
        if (!hidden && model) {
            model->draw();
        }
    }
    void render_line_model() {
        if (!hidden && model) {
            directionModel->draw();
        }
    }

    void render_properties(){

        ImGui::PushID(this);

        bool changedRadius = false;
        bool changedMetric = false;
        bool changedDirection = false;

        ImGui::ColorEdit4("Color", (float*)&color);
        ImGui::ColorEdit4("Line Color", (float*)&lineColor);

        changedDirection |= ImGui::InputFloat3("Origin", origin);
        
        changedRadius |= ImGui::InputFloat(
            "Sigma (influence radius ≈ 3×sigma)",
            &sigma, 0.1f, 50.0f, "%.2f");

        changedMetric |= ImGui::InputFloat(
            "Stretch X (along Dir)",
            &stretch.x,
            0.01f, 5.0f, "%.3f"
        );
        changedMetric |=ImGui::InputFloat(
            "Stretch Y (transverse)",
            &stretch.y,
            0.01f, 5.0f, "%.3f"
        );
        changedMetric |=ImGui::InputFloat(
            "Stretch Z (transverse)",
            &stretch.z,
            0.01f, 5.0f, "%.3f"
        );
        
        changedDirection |= ImGui::InputFloat(
            "Roll Angle", &angle, 0.01f, 10.0f, "%.4f");
        
            ImGui::InputFloat3(
            "Material Direction (local X)", direction);
        
        ImGui::SameLine();
        if (ImGui::Button("Normalize")){
            if (direction.norm() > 1e-5f){
                direction.normalize();
                changedDirection = true;
            }
        }
        
        if (changedMetric || changedDirection){
            update_metric();
        }

        if (changedRadius || changedDirection){
            update_model();
        }

        ImGui::PopID();

    };

    ObjectType get_type() const {
		return ObjectType::AnisoSource;
	}

    void update_metric() {

        Eigen::Matrix3f rot = rotation_from_direction_roll(
            direction, 
            angle
        );
        
        auto safe_sq = [](float val) { return std::max(1e-6f, val * val); };
        Eigen::Vector3f invSq(1.0f / safe_sq(stretch.x),
                              1.0f / safe_sq(stretch.y),
                              1.0f / safe_sq(stretch.z));
        M = rot.transpose() * invSq.asDiagonal() * rot;

    };

    void update_model(){
        float radius = 3.0f * sigma;
        Vec3 end = radius * direction;
        if (!model) {
            model = std::make_unique<Sphere>(radius, 36, 18);
            directionModel = std::make_unique<LineModel>(
                Vec3{0.0f, 0.0f, 0.f}, end);
        } else {
            model->update_radius(radius);
            directionModel = std::make_unique<LineModel>(
                Vec3{0.0f, 0.0f, 0.f}, end);
        }
    };

private:
    std::unique_ptr<Sphere> model;
    std::unique_ptr<LineModel> directionModel;
};

//@brief function to create the tensor M for a single source
inline void create_metric(AnisotropySource& source){
    // rotation matrix is using the applied direction
    Eigen::Matrix3f rot = rotation_from_direction_roll(
        source.direction, source.angle
    );
    
    // S^-2 matrix diagonal
    Eigen::Vector3f invSq(1.0f / (source.stretch.x * source.stretch.x),
                          1.0f / (source.stretch.y * source.stretch.y),
                          1.0f / (source.stretch.z * source.stretch.z));
    // M = R^T S^-2 R
    source.M = rot.transpose() * invSq.asDiagonal() * rot;
};

//@ brief function to estimate the number of candidates based on the applied stretches
inline size_t choose_candidate_number(
    float stretchMin, float stretchMax, size_t seedCount)
{
    if (seedCount == 0) return 0;
    float  kappa = stretchMax / std::max(stretchMin, 1e-6f);
    size_t k     = static_cast<size_t>(std::ceil(3.0f * kappa * kappa));
    size_t lo    = std::min<size_t>(12, seedCount);
    k = std::max(k, lo);
    k = std::min(k, seedCount);
    return k;
}

inline Eigen::Matrix3f blend_metric(
    const Vec3 point, 
    const std::vector<std::shared_ptr<AnisotropySource>>& sources,
    AnisotropySource& background, 
    float backgroundWeight){

        Eigen::Matrix3f M = backgroundWeight * background.M;
        float weightSum = backgroundWeight;
        for(const auto& src: sources){
            Vec3 d = point - src->origin;
            float r2 = d.x * d.x + d.y * d.y + d.z * d.z;
            if (r2 > 9.0f * src->sigma * src->sigma) continue; 

            float invSigma = 1.0f / (2.0f * src->sigma * src->sigma);
            float w = std::exp(-r2 * invSigma);
            M += w * src->M;
            weightSum += w;
        }

        // return the weighted sum
        return M / weightSum;
    };

inline double aniso_distance_sq(
    const Eigen::Matrix3f& M,
    const Vec3& p, 
    const Vec3& q){
    Eigen::Vector3f delta(p.x - q.x, p.y - q.y, p.z - q.z);
    return static_cast<double>(delta.dot(M * delta));
};

inline Vec3 aniso_distance_grad(
    const Eigen::Matrix3f& M,
    const Vec3& p, const Vec3& q, float d)
{
    if (d < 1e-6f) return Vec3(0.0f, 0.0f, 0.0f);  // at the seed centre
    Eigen::Vector3f delta(p.x - q.x, p.y - q.y, p.z - q.z);
    Eigen::Vector3f g = (M * delta) / d;
    return Vec3(g.x(), g.y(), g.z());
}
#endif