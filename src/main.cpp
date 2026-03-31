#include "../include/MeshData.hpp"
#include "../include/AstroProperties.hpp"
#include "../include/Boulder_detector.hpp"

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
//             MeshData data = getMeshData(mesh, density);
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
        boulder::PipelineParams params;
        params.mesh_path = "/Users/thangthend11/Desktop/Hackathon-GIG-2026/data/bennu_OLA_v21_PTM_very-high.obj";
        params.output_dir = "boulder_demo_single_scale";

        params.smooth_iterations = 0;
        params.radius_multiplier = 6.0;
        params.min_neighbors = 20;
        params.use_absolute_height = false;

        params.height_threshold = std::numeric_limits<double>::quiet_NaN();
        params.height_sigma_factor = 1.5;
        params.roughness_max = std::numeric_limits<double>::quiet_NaN();

        params.outer_ring_depth = 1;

        params.min_vertices = 30;
        params.min_max_height = 0.0;
        params.min_diameter = 0.0;
        params.max_diameter = std::numeric_limits<double>::infinity();
        params.min_compactness = 0.0;
        params.max_elongation = 3.0;
        params.min_prominence_mean = 0.0;
        params.min_prominence_max = 0.0;
        params.min_boundary_drop = 0.0;
        params.min_boundary_irregularity = 0.0;
        params.max_boundary_irregularity = std::numeric_limits<double>::infinity();

        params.visualize = true;
        params.verbose = true;

        auto result = boulder::DetectBouldersSingleScale(params);
        std::cout << "[INFO] Final kept boulders = " << result.kept_boulders.size() << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
    return 0;
}