#pragma once

#include <open3d/Open3D.h>
#include <Eigen/Dense>

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

namespace boulder {

struct PlaneFitResult {
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    Eigen::Vector3d evals = Eigen::Vector3d::Zero();  // ascending
    Eigen::Matrix3d evecs = Eigen::Matrix3d::Identity();
};

struct LocalHeightResult {
    std::vector<double> heights;
    std::vector<std::uint8_t> valid_mask;
    std::vector<double> local_roughness;
    std::vector<double> local_planarity;
    std::vector<int> neighbor_counts;
};

struct CandidateResult {
    std::vector<std::uint8_t> candidate_mask;
    double used_threshold = 0.0;
    double mu = 0.0;
    double sigma = 0.0;
};

struct BoundaryMetrics {
    std::vector<int> boundary_vertices;
    std::vector<int> outer_ring_vertices;

    double inside_mean_height = 0.0;
    double inside_max_height = 0.0;
    double boundary_mean_height = 0.0;
    double boundary_height_std = 0.0;
    double outer_ring_mean_height = 0.0;
    double outer_ring_height_std = 0.0;

    double prominence_mean = 0.0;
    double prominence_max = 0.0;
    double boundary_drop = 0.0;
    double boundary_irregularity = 0.0;

    int boundary_size = 0;
    int outer_ring_size = 0;
};

struct ComponentInfo {
    int component_id = -1;
    std::vector<int> vertex_indices;
    int num_vertices = 0;
    std::array<double, 3> centroid = {0.0, 0.0, 0.0};

    double max_height = 0.0;
    double mean_height = 0.0;
    double diameter_estimate = 0.0;
    double horizontal_size = 0.0;
    double surface_area_estimate = 0.0;
    double compactness_ratio = 0.0;
    double elongation = 1.0;
    std::vector<double> pca_eigenvalues;

    double inside_mean_height = 0.0;
    double inside_max_height = 0.0;
    double boundary_mean_height = 0.0;
    double boundary_height_std = 0.0;
    double outer_ring_mean_height = 0.0;
    double outer_ring_height_std = 0.0;

    double prominence_mean = 0.0;
    double prominence_max = 0.0;
    double boundary_drop = 0.0;
    double boundary_irregularity = 0.0;

    int boundary_size = 0;
    int outer_ring_size = 0;
    std::vector<int> boundary_vertices;
    std::vector<int> outer_ring_vertices;
};

struct DetectionResult {
    std::shared_ptr<open3d::geometry::TriangleMesh> mesh;
    double mean_edge_length = 0.0;
    double radius = 0.0;

    std::vector<double> heights;
    std::vector<std::uint8_t> valid_mask;
    std::vector<std::uint8_t> candidate_mask;
    std::vector<std::uint8_t> pred_vertex_mask;
    std::vector<int> labels;
    std::vector<std::vector<int>> components;

    std::vector<ComponentInfo> all_component_infos;
    std::vector<ComponentInfo> kept_boulders;
    std::vector<ComponentInfo> rejected_components;

    double used_threshold = 0.0;
};

struct BoulderSizeDistribution {
        std::vector<double> diameters_sorted;
        std::vector<double> bin_edges;
        std::vector<double> bin_centers;
        std::vector<int> bin_counts;
        std::vector<int> cumulative_counts_ge;
        std::vector<double> cumulative_number_density_ge;
};

struct PipelineParams {
    std::string mesh_path;
    std::string output_dir = "boulder_output_single_scale";

    int smooth_iterations = 0;
    double radius_multiplier = 6.0;
    int min_neighbors = 20;
    bool use_absolute_height = false;

    double height_threshold = std::numeric_limits<double>::quiet_NaN();
    double height_sigma_factor = 1.5;
    double roughness_max = std::numeric_limits<double>::quiet_NaN();
    int outer_ring_depth = 1;

    int min_vertices = 30;
    double min_max_height = 0.0;
    double min_diameter = 0.0;
    double max_diameter = std::numeric_limits<double>::infinity();
    double min_compactness = 0.0;
    double max_elongation = 3.0;
    double min_prominence_mean = 0.0;
    double min_prominence_max = 0.0;
    double min_boundary_drop = 0.0;
    double min_boundary_irregularity = 0.0;
    double max_boundary_irregularity = std::numeric_limits<double>::infinity();

    bool visualize = true;
    bool verbose = true;

    //for merge 
    bool enable_component_merging = false;
    int merge_knn = 24;
    double merge_boundary_distance_factor = 4.0;
    double merge_centroid_distance_factor = 1.5;
    double merge_height_difference_factor = 2.0;
    double merge_prominence_difference_factor = 2.0;
};


std::shared_ptr<open3d::geometry::TriangleMesh> LoadMesh(const std::string& mesh_path,
                                                         bool verbose = true);

std::shared_ptr<open3d::geometry::TriangleMesh> OptionalSmoothMesh(
        const open3d::geometry::TriangleMesh& mesh,
        int iterations = 0);

double ComputeMeanEdgeLength(const open3d::geometry::TriangleMesh& mesh);

std::vector<std::vector<int>> BuildTopologyAdjacency(
        const open3d::geometry::TriangleMesh& mesh);

std::vector<double> TriangleAreas(const open3d::geometry::TriangleMesh& mesh);

std::vector<std::vector<int>> BuildVertexToTriangleMap(
        const open3d::geometry::TriangleMesh& mesh);

PlaneFitResult FitPlanePCA(const std::vector<Eigen::Vector3d>& points);

LocalHeightResult ComputeLocalHeights(const open3d::geometry::TriangleMesh& mesh,
                                      double radius,
                                      int min_neighbors = 20,
                                      bool use_absolute_height = false,
                                      bool verbose = true);

CandidateResult DetectCandidateVertices(const std::vector<double>& heights,
                                        const std::vector<std::uint8_t>& valid_mask,
                                        const std::vector<double>* local_roughness = nullptr,
                                        double height_threshold = std::numeric_limits<double>::quiet_NaN(),
                                        double height_sigma_factor = 1.5,
                                        double roughness_max = std::numeric_limits<double>::quiet_NaN(),
                                        bool verbose = true);

void ExtractCandidateComponents(const std::vector<std::uint8_t>& candidate_mask,
                                const std::vector<std::vector<int>>& adjacency,
                                std::vector<int>& labels,
                                std::vector<std::vector<int>>& components,
                                bool verbose = true);

double EstimateComponentDiameter(const std::vector<Eigen::Vector3d>& points);

double EstimateComponentHorizontalSize(const std::vector<Eigen::Vector3d>& points);

double ComputeComponentArea(const std::vector<int>& component_vertices,
                            const std::vector<std::vector<int>>& vertex_to_triangles,
                            const std::vector<double>& tri_areas);

std::pair<double, std::vector<double>> ComputeComponentPCAShape(
        const std::vector<Eigen::Vector3d>& points);

void GetComponentBoundaryAndOuterRing(const std::vector<int>& component_vertices,
                                      const std::vector<std::vector<int>>& adjacency,
                                      int outer_ring_depth,
                                      std::vector<int>& boundary,
                                      std::vector<int>& outer_ring);

BoundaryMetrics ComputeBoundaryMetrics(const std::vector<int>& component_vertices,
                                       const std::vector<double>& heights,
                                       const std::vector<std::vector<int>>& adjacency,
                                       int outer_ring_depth = 1);

std::vector<ComponentInfo> MeasureComponents(
        const open3d::geometry::TriangleMesh& mesh,
        const std::vector<std::vector<int>>& components,
        const std::vector<double>& heights,
        const std::vector<std::vector<int>>& adjacency,
        const std::vector<std::vector<int>>& vertex_to_triangles,
        const std::vector<double>& tri_areas,
        int outer_ring_depth = 1,
        bool verbose = true);

void FilterBoulderComponents(const std::vector<ComponentInfo>& component_infos,
                             std::vector<ComponentInfo>& kept,
                             std::vector<ComponentInfo>& rejected,
                             int min_vertices = 30,
                             double min_max_height = 0.0,
                             double min_diameter = 0.0,
                             double max_diameter = std::numeric_limits<double>::infinity(),
                             double min_compactness = 0.0,
                             double max_elongation = 3.0,
                             double min_prominence_mean = 0.0,
                             double min_prominence_max = 0.0,
                             double min_boundary_drop = 0.0,
                             double min_boundary_irregularity = 0.0,
                             double max_boundary_irregularity = std::numeric_limits<double>::infinity(),
                             bool verbose = true);

BoulderSizeDistribution ComputeBoulderSizeDistribution(
                        const std::vector<ComponentInfo>& kept_boulders,
                        double total_surface_area,
                        int num_bins = 20,
                        bool log_bins = true,
                        bool verbose = true);


std::vector<Eigen::Vector3d> ScalarToColors(const std::vector<double>& values,
                                            const std::vector<std::uint8_t>* valid_mask = nullptr);

std::shared_ptr<open3d::geometry::TriangleMesh> ColorizeScalarField(
        const open3d::geometry::TriangleMesh& mesh,
        const std::vector<double>& scalars,
        const std::vector<std::uint8_t>* valid_mask = nullptr);

std::shared_ptr<open3d::geometry::TriangleMesh> ColorizeCandidateVertices(
        const open3d::geometry::TriangleMesh& mesh,
        const std::vector<std::uint8_t>& candidate_mask,
        const Eigen::Vector3d& base_color = Eigen::Vector3d(0.7, 0.7, 0.7),
        const Eigen::Vector3d& candidate_color = Eigen::Vector3d(1.0, 0.0, 0.0));

std::shared_ptr<open3d::geometry::TriangleMesh> ColorizeComponents(
        const open3d::geometry::TriangleMesh& mesh,
        const std::vector<ComponentInfo>& component_infos);

bool SaveMesh(const open3d::geometry::TriangleMesh& mesh, const std::string& output_path);
void SaveComponentInfosJson(const std::vector<ComponentInfo>& component_infos,
                            const std::string& output_json_path);
void SaveDoubleVectorTxt(const std::vector<double>& values, const std::string& path);
void SaveIntVectorTxt(const std::vector<int>& values, const std::string& path);
void SaveMaskTxt(const std::vector<std::uint8_t>& values, const std::string& path);

DetectionResult DetectBouldersSingleScale(const PipelineParams& params);
std::vector<std::vector<int>> BuildComponentMergeGroups(
        const open3d::geometry::TriangleMesh& mesh,
        const std::vector<ComponentInfo>& component_infos,
        double mean_edge_length,
        int merge_knn = 24,
        double merge_boundary_distance_factor = 4.0,
        double merge_centroid_distance_factor = 1.5,
        double merge_height_difference_factor = 2.0,
        double merge_prominence_difference_factor = 2.0,
        bool verbose = true);

std::vector<std::vector<int>> BuildMergedComponents(
        const std::vector<ComponentInfo>& component_infos,
        const std::vector<std::vector<int>>& merge_groups);

}  // namespace boulder

