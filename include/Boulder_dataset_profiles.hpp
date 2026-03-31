#pragma once
#include "Boulder_detector.hpp"
#include "Boulder_profiles.hpp"
#include <limits>

namespace dataset_profiles {

inline boulder::PipelineParams RyuguPersonalHiRes(const std::string& mesh_path,
                                                  const std::string& output_dir) {
    auto p = boulder_profiles::Ryugu(mesh_path, output_dir);

    auto mesh = boulder::LoadMesh(mesh_path, p.verbose);
    const double e = boulder::ComputeMeanEdgeLength(*mesh);

    p.smooth_iterations = 0;
    p.radius_multiplier = 5.8;
    p.min_neighbors = 20;

    p.height_threshold = std::numeric_limits<double>::quiet_NaN();
    p.height_sigma_factor = 1.35;
    p.roughness_max = std::numeric_limits<double>::quiet_NaN();
    p.outer_ring_depth = 1;

    p.min_vertices = 24;
    p.min_max_height = 0.8 * e;
    p.min_diameter = 3.5 * e;
    p.max_diameter = std::numeric_limits<double>::infinity();

    p.min_compactness = 0.005;
    p.max_elongation = 3.2;

    p.min_prominence_mean = 0.25 * e;
    p.min_prominence_max  = 0.80 * e;
    p.min_boundary_drop   = 0.20 * e;

    p.min_boundary_irregularity = 0.0;
    p.max_boundary_irregularity = std::numeric_limits<double>::infinity();

    p.enable_component_merging = true;
    p.merge_knn = 20;
    p.merge_boundary_distance_factor = 3.5;
    p.merge_centroid_distance_factor = 1.4;
    p.merge_height_difference_factor = 1.8;
    p.merge_prominence_difference_factor = 1.8;

    return p;
}
inline boulder::PipelineParams BennuPersonalHiRes(const std::string& mesh_path,
                                                  const std::string& output_dir) {
    auto p = boulder_profiles::Bennu(mesh_path, output_dir);

    auto mesh = boulder::LoadMesh(mesh_path, p.verbose);
    const double e = boulder::ComputeMeanEdgeLength(*mesh);

    p.min_vertices = 22;

    p.min_max_height = 0.7 * e;
    p.min_diameter   = 3.0 * e;
    p.max_diameter   = std::numeric_limits<double>::infinity();

    p.min_compactness = 0.003;

    p.min_prominence_mean = 0.20 * e;
    p.min_prominence_max  = 0.70 * e;
    p.min_boundary_drop   = 0.16 * e;

    p.min_boundary_irregularity = 0.0;
    p.max_boundary_irregularity = std::numeric_limits<double>::infinity();

    p.enable_component_merging = true;
    p.merge_knn = 20;
    p.merge_boundary_distance_factor = 3.8;
    p.merge_centroid_distance_factor = 1.5;
    p.merge_height_difference_factor = 1.8;
    p.merge_prominence_difference_factor = 1.8;

    return p;
}

inline boulder::PipelineParams ChuryPersonalHiRes(const std::string& mesh_path,
                                                  const std::string& output_dir) {
    auto p = boulder_profiles::Chury(mesh_path, output_dir);

    auto mesh = boulder::LoadMesh(mesh_path, p.verbose);
    const double e = boulder::ComputeMeanEdgeLength(*mesh);

    p.min_vertices = 32;

    p.min_max_height = 1.2 * e;
    p.min_diameter   = 4.0 * e;
    p.max_diameter   = std::numeric_limits<double>::infinity();

    p.min_compactness = 0.008;

    p.min_prominence_mean = 0.22 * e;
    p.min_prominence_max  = 0.75 * e;
    p.min_boundary_drop   = 0.18 * e;

    p.min_boundary_irregularity = 0.0;
    p.max_boundary_irregularity = std::numeric_limits<double>::infinity();

    p.enable_component_merging = true;
    p.merge_knn = 20;
    p.merge_boundary_distance_factor = 3.6;
    p.merge_centroid_distance_factor = 1.35;
    p.merge_height_difference_factor = 1.7;
    p.merge_prominence_difference_factor = 1.7;

    return p;
}

inline boulder::PipelineParams PhobosPersonalHiRes(const std::string& mesh_path,
                                                   const std::string& output_dir) {
    auto p = boulder_profiles::Phobos(mesh_path, output_dir);

    auto mesh = boulder::LoadMesh(mesh_path, p.verbose);
    const double e = boulder::ComputeMeanEdgeLength(*mesh);

    p.min_vertices = 28;

    p.min_max_height = 1.0 * e;
    p.min_diameter   = 3.8 * e;
    p.max_diameter   = std::numeric_limits<double>::infinity();

    p.min_compactness = 0.010;

    p.min_prominence_mean = 0.30 * e;
    p.min_prominence_max  = 0.90 * e;
    p.min_boundary_drop   = 0.25 * e;

    p.min_boundary_irregularity = 0.0;
    p.max_boundary_irregularity = std::numeric_limits<double>::infinity();

    p.enable_component_merging = true;
    p.merge_knn = 18;
    p.merge_boundary_distance_factor = 3.2;
    p.merge_centroid_distance_factor = 1.3;
    p.merge_height_difference_factor = 1.6;
    p.merge_prominence_difference_factor = 1.6;

    return p;
}

}