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

    struct InputTextCallback_UserData
    {
        std::string* Str;
        ImGuiInputTextCallback  ChainCallback;
        void* ChainCallbackUserData;
    };

    static int InputTextCallback(ImGuiInputTextCallbackData* data)
    {
        InputTextCallback_UserData* user_data = (InputTextCallback_UserData*)data->UserData;
        if (data->EventFlag == ImGuiInputTextFlags_CallbackResize)
        {
            // Resize string callback
            // If for some reason we refuse the new length (BufTextLen) and/or capacity (BufSize) we need to set them back to what we want.
            std::string* str = user_data->Str;
            IM_ASSERT(data->Buf == str->c_str());
            str->resize(data->BufTextLen);
            data->Buf = (char*)str->c_str();
        }
        else if (user_data->ChainCallback)
        {
            // Forward to user callback, if any
            data->UserData = user_data->ChainCallbackUserData;
            return user_data->ChainCallback(data);
        }
        return 0;
    }

    bool ImGui::InputText(const char* label, std::string* str, ImGuiInputTextFlags flags, ImGuiInputTextCallback callback, void* user_data)
    {
        IM_ASSERT((flags & ImGuiInputTextFlags_CallbackResize) == 0);
        flags |= ImGuiInputTextFlags_CallbackResize;

        InputTextCallback_UserData cb_user_data;
        cb_user_data.Str = str;
        cb_user_data.ChainCallback = callback;
        cb_user_data.ChainCallbackUserData = user_data;
        return InputText(label, (char*)str->c_str(), str->capacity() + 1, flags, InputTextCallback, &cb_user_data);
    }

    bool ImGui::InputTextMultiline(const char* label, std::string* str, const ImVec2& size, ImGuiInputTextFlags flags, ImGuiInputTextCallback callback, void* user_data)
    {
        IM_ASSERT((flags & ImGuiInputTextFlags_CallbackResize) == 0);
        flags |= ImGuiInputTextFlags_CallbackResize;

        InputTextCallback_UserData cb_user_data;
        cb_user_data.Str = str;
        cb_user_data.ChainCallback = callback;
        cb_user_data.ChainCallbackUserData = user_data;
        return InputTextMultiline(label, (char*)str->c_str(), str->capacity() + 1, size, flags, InputTextCallback, &cb_user_data);
    }

    bool ImGui::InputTextWithHint(const char* label, const char* hint, std::string* str, ImGuiInputTextFlags flags, ImGuiInputTextCallback callback, void* user_data)
    {
        IM_ASSERT((flags & ImGuiInputTextFlags_CallbackResize) == 0);
        flags |= ImGuiInputTextFlags_CallbackResize;

        InputTextCallback_UserData cb_user_data;
        cb_user_data.Str = str;
        cb_user_data.ChainCallback = callback;
        cb_user_data.ChainCallbackUserData = user_data;
        return InputTextWithHint(label, hint, (char*)str->c_str(), str->capacity() + 1, flags, InputTextCallback, &cb_user_data);
    }

}
