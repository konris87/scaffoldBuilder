#include "Anisotropy.h"

void AnisotropyFactory::launch(){

    name = "";
    pendingSrc = std::make_shared<AnisotropySource>();
};

void AnisotropyFactory::render(
    Logger* logger,
    const char* popupName, bool& showPopup,
    std::vector<std::shared_ptr<AnisotropySource>>& globalSources
){

    if (showPopup) {
		ImGui::OpenPopup(popupName);
	}

	// always centered
	ImVec2 center = ImGui::GetMainViewport()->GetCenter();
	ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

	if (ImGui::BeginPopupModal(popupName, &showPopup, ImGuiWindowFlags_AlwaysAutoResize))
	{
		// set the source name
		ImGui::InputText("Name", &name);

		ImGui::SeparatorText("Properties");
        pendingSrc->render_properties();

		if(ImGui::Button("Create")){
			pendingSrc->update_metric();
			pendingSrc->update_model();
			if (name.empty()){
				pendingSrc->name = "Anisotropy Source " + std::to_string(globalSources.size() + 1);
			}
            else{
                pendingSrc->name = name;
            }
			globalSources.push_back(std::move(pendingSrc));
            showPopup = false;
            ImGui::EndPopup();
		}
		ImGui::SameLine();
		if (ImGui::Button("Cancel")){
            showPopup = false;
            ImGui::EndPopup();
		}

		ImGui::EndPopup();
	}

};