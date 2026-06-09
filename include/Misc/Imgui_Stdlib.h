#ifndef IMGUI_STDLIB_H
#define IMGUI_STDLIB_H

#include <array>
#include <Math/Vec.h>

namespace ImGui
{
    // ImGui::InputFloat3() with std::array
    IMGUI_API bool InputFloat3(const char* label, std::array<float, 3>& val, const char* format = "%.3f", ImGuiInputTextFlags flags = 0);

    IMGUI_API bool InputFloat3(
        const char* label, Vec3& val, const char* format = "%.3f", ImGuiInputTextFlags flags = 0);

    IMGUI_API bool InputFloat2(const char* label, std::array<float, 2>& val, const char* format = "%.3f", ImGuiInputTextFlags flags = 0);

    IMGUI_API bool SliderFloat3(
        const char* label, std::array<float, 3>& val, float v_min, float v_max, const char* format = "%.3f", ImGuiInputTextFlags flags = 0);

    IMGUI_API bool InputInt3(const char* label, std::array<int, 3>& val, ImGuiInputTextFlags flags = 0);

    IMGUI_API bool DragFloat3(const char* label, Vec3& val, float v_speed, float v_min, float v_max, const char* format = "%.3f", ImGuiInputTextFlags flags = 0);

	IMGUI_API bool SetDragDropPayload(const char* type, const Vec3& data, ImGuiCond cond = 0);

    IMGUI_API bool  InputText(const char* label, std::string* str, ImGuiInputTextFlags flags = 0, ImGuiInputTextCallback callback = nullptr, void* user_data = nullptr);

    IMGUI_API bool  InputTextMultiline(const char* label, std::string* str, const ImVec2& size = ImVec2(0, 0), ImGuiInputTextFlags flags = 0, ImGuiInputTextCallback callback = nullptr, void* user_data = nullptr);

    IMGUI_API bool  InputTextWithHint(const char* label, const char* hint, std::string* str, ImGuiInputTextFlags flags = 0, ImGuiInputTextCallback callback = nullptr, void* user_data = nullptr);
}

#endif