#include "imgui.h"
#include "Imgui_Stdlib.h"
namespace ImGui
{

    bool ImGui::InputFloat3(const char* label, std::array<float, 3>& val, const char* format, ImGuiInputTextFlags flags) {

        return InputScalarN(label, ImGuiDataType_Float, val.data(), 3, NULL, NULL, format, flags);
    };

    bool ImGui::InputFloat2(const char* label, std::array<float, 2>& val, const char* format, ImGuiInputTextFlags flags) {

        return InputScalarN(label, ImGuiDataType_Float, val.data(), 2, NULL, NULL, format, flags);
    };


    bool ImGui::SliderFloat3(const char* label, std::array<float, 3>& val, float v_min, float v_max, const char* format, ImGuiInputTextFlags flags) {

        return SliderScalarN(label, ImGuiDataType_Float, val.data(), 3, &v_min, &v_max, format, flags);
    };


    bool ImGui::InputInt3(const char* label, std::array<int, 3>& val, ImGuiInputTextFlags flags) {

        return InputScalarN(label, ImGuiDataType_S32, val.data(), 3, NULL, NULL, "%d", flags);
    };


    // overload Vec class
    bool InputFloat3(const char* label, Vec3& v,
        const char* format, ImGuiInputTextFlags flags)
    {
        return InputScalarN(label,
            ImGuiDataType_Float,
            &v.x,  // pointer to first component
            3,
            nullptr, nullptr,
            format,
            flags);
    }

    // overload dragable Vec3
        // overload Vec class
    bool DragFloat3(
        const char* label, Vec3& val,
        float v_speed, float v_min,
        float v_max, const char* format, ImGuiInputTextFlags flags)
    {
        return DragScalarN(
            label,
            ImGuiDataType_Float,
            &val.x,  // pointer to first component
            3,
            v_speed, &v_min, &v_max,
            format,
            flags);
    }

    bool SetDragDropPayload(const char* type, const Vec3& data, ImGuiCond cond) {
    
		return SetDragDropPayload(type, &data.x, sizeof(float) * 3, cond);
    };

}
