#include "../include/MeshData.hpp"
#include "../include/AstroProperties.hpp"
#include "../include/Boulder_detector.hpp"
#include "../include/Boulder_profiles.hpp"
#include "../include/Boulder_dataset_profiles.hpp"
#include <iostream>
#include <string>
#include <vector>

struct AstroInput {
    std::string name;
    std::string mesh_path;
};

// int main() {
//     // Sửa các path này theo thư mục dữ liệu thực tế của bạn
//     std::vector<AstroInput> astros = {
//         {"Chury",   "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/chury-dlr_spg-shap7-v1.0.ply"},
//         {"Phobos",  "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/g_phobos_018m_spc_0000n00000_v002.obj"},
//         {"Bennu",   "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/bennu_OLA_v21_PTM_very-high.obj"},
//         {"Ryugu",   "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/SHAPE_SFM_RYUGU_3M_v20180804.obj"},
//         {"Lune",    "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/moon_ldem_64-cart-020m.ply"},
//         {"Mars",    "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/mars-megr90n000fb-cart-020m.ply"},
//         {"Mercure", "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/MERCURY_MSGR_DEM_USG_SC_I_V02-cart-020m.ply"}
//     };

//     for (const auto& astro : astros) {
//         try {
//             Mesh mesh;
//             if (!loadMesh(astro.mesh_path, mesh)) {
//                 std::cerr << "Skipping " << astro.name << '\n';
//                 continue;
//             }

//             double density = getDensityKgPerM3(astro.name);
//             MeshData data = getMeshData(mesh, density, LinearUnit::Kilometer);
//             displayMeshData(astro.name, data);
//         }
//         catch (const std::exception& e) {
//             std::cerr << "Error for " << astro.name << ": " << e.what() << '\n';
//         }
//     }

//     return 0;
// }

int main() {
    try {
        auto params = dataset_profiles::RyuguPersonalHiRes(
            "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/SHAPE_SFM_RYUGU_3M_v20180804.obj",
            "boulder_demo_ryugu"
        );

        // Ví dụ đổi sang:
        // auto params = dataset_profiles::RyuguPersonalHiRes(...);
        // auto params = dataset_profiles::ChuryPersonalHiRes(...);
        // auto params = dataset_profiles::PhobosPersonalHiRes(...);

        auto result = boulder::DetectBouldersSingleScale(params);

        std::cout << "[INFO] Final kept boulders = "
                  << result.kept_boulders.size() << "\n";
        std::cout << "[INFO] mean_edge_length    = "
                  << result.mean_edge_length << "\n";
        std::cout << "[INFO] radius used        = "
                  << result.radius << "\n";
        std::cout << "[INFO] threshold used     = "
                  << result.used_threshold << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
    return 0;
}