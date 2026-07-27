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
    ) : origin(center), direction(direction), angle(angle), stretch(stretch), sigma(sigma){
        update_metric();
    };

    Vec3 origin; // set initially to zero
    Vec3 direction = {1.0f, 0.0f, 0.0f};
    float angle = 0.0f;
    Vec3 stretch = {1.0f, 1.0f, 1.0f};
    float sigma = 1.1f; // gaussian falloff parameter
    Eigen::Matrix3f C = Eigen::Matrix3f::Identity();

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

        // Use Covariance Matrix
        Eigen::Matrix3f rot = rotation_from_direction_roll(
            direction, 
            angle
        );
        
        auto safe_sq = [](float val) { return std::max(1e-6f, val * val); };

        Eigen::Vector3f sq(
            safe_sq(stretch.x), 
            safe_sq(stretch.y), 
            safe_sq(stretch.z)
        );

        C = rot.transpose() * sq.asDiagonal() * rot;
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

    auto safe_sq = [](float val) { return std::max(1e-6f, val * val); };
    
    Eigen::Vector3f sq(
        safe_sq(source.stretch.x),
        safe_sq(source.stretch.y),
        safe_sq(source.stretch.z)
    );

    // Covariance matrix R^T S^2 R
    source.C = rot.transpose() * sq.asDiagonal() * rot;
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

//@brief Blend the source covariances at 'point' with the global background
// covariance, then invert to obtain the local Riemannian metric M(x).
//
// Partition of unity: the background acts as an always-present source with
// constant weight wb, so the blend is a single normalized average
//     C(x) = ( sum_i w_i C_i + wb * C_bg ) / ( sum_i w_i + wb ).
// Far from every source C(x) = C_bg exactly; at a source centre the source
// contributes 1/(1+wb) of the blend. Larger wb bleeds the background more
// strongly into source regions (wb = 1 caps a lone source at a 50/50 mix);
// wb -> 0 would make the blend undefined outside all supports, so it is
// clamped away from zero.
//
// The Gaussian weight is shifted so it reaches zero continuously at the
// 3-sigma support boundary (an unshifted kernel jumps by exp(-4.5) ~ 0.011
// there, leaving a small discontinuity in the scalar field).
inline Eigen::Matrix3f blend_metric(
    const Vec3 point,
    const std::vector<std::shared_ptr<AnisotropySource>>& sources,
    const AnisotropySource& background,
    float backgroundWeight = 0.1f){

        // weight of the unshifted Gaussian at the 3-sigma cutoff
        const float cutoffWeight = std::exp(-4.5f);

        Eigen::Matrix3f Clocal = Eigen::Matrix3f::Zero();
        float weightSum = 0.0f;
        for(const auto& src : sources){

            // squared distance of query point from center of source
            Vec3 d = point - src->origin;

            float r2 = d.x * d.x + d.y * d.y + d.z * d.z;
            if (r2 > 9.0f * src->sigma * src->sigma) continue;

            // Gaussian falloff, shifted to vanish at the support boundary
            float invSigma = 1.0f / (2.0f * src->sigma * src->sigma);
            float w = std::exp(-r2 * invSigma) - cutoffWeight;
            if (w <= 0.0f) continue;

            Clocal += w * src->C;
            weightSum += w;
        }

        float wb = std::max(backgroundWeight, 1e-3f);
        Eigen::Matrix3f Cblend = (Clocal + wb * background.C) / (weightSum + wb);
        return Cblend.inverse();
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