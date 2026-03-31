#include "../include/Boulder_detector.hpp"

//#include <Open3D/IO/ClassIO/TriangleMeshIO.h>
#include <Open3D/Visualization/Utility/DrawGeometry.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

namespace boulder {
namespace {

using open3d::geometry::KDTreeFlann;
using open3d::geometry::PointCloud;
using open3d::geometry::TriangleMesh;

constexpr double kEps = 1e-12;

struct DisjointSet {
    std::vector<int> parent;
    std::vector<int> rank;

    explicit DisjointSet(int n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int Find(int x) {
        if (parent[x] != x) {
            parent[x] = Find(parent[x]);
        }
        return parent[x];
    }

    void Union(int a, int b) {
        a = Find(a);
        b = Find(b);
        if (a == b) return;

        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
    }
};

struct CsvRowInt {
    double x = 0.0;
    int y = 0;
};

struct CsvRowDouble {
    double x = 0.0;
    double y = 0.0;
};

void SaveCsvXY_Int(const std::vector<double>& xs,
                   const std::vector<int>& ys,
                   const std::string& path,
                   const std::string& header_x,
                   const std::string& header_y) {
    std::ofstream ofs(path);
    ofs << header_x << "," << header_y << "\n";
    const std::size_t n = std::min(xs.size(), ys.size());
    ofs << std::setprecision(16);
    for (std::size_t i = 0; i < n; ++i) {
        ofs << xs[i] << "," << ys[i] << "\n";
    }
}

void SaveCsvXY_Double(const std::vector<double>& xs,
                      const std::vector<double>& ys,
                      const std::string& path,
                      const std::string& header_x,
                      const std::string& header_y) {
    std::ofstream ofs(path);
    ofs << header_x << "," << header_y << "\n";
    const std::size_t n = std::min(xs.size(), ys.size());
    ofs << std::setprecision(16);
    for (std::size_t i = 0; i < n; ++i) {
        ofs << xs[i] << "," << ys[i] << "\n";
    }
}

void WriteBoulderDistributionGnuplotScript(const std::string& output_dir,
                                           const std::string& body_label) {
    std::ofstream gp(output_dir + "/plot_boulder_distribution.gp");

    gp << "set datafile separator ','\n";
    gp << "set terminal pngcairo size 1400,700\n";

    gp << "set output '" << output_dir << "/boulder_size_histogram.png'\n";
    gp << "set title 'Boulder size histogram - " << body_label << "'\n";
    gp << "set xlabel 'Estimated diameter'\n";
    gp << "set ylabel 'Count'\n";
    gp << "set grid\n";
    gp << "plot '" << output_dir
       << "/boulder_size_histogram.csv' using 1:2 with boxes title 'Histogram'\n";

    gp << "set output '" << output_dir << "/boulder_size_cumulative.png'\n";
    gp << "set title 'Boulder cumulative size-frequency - " << body_label << "'\n";
    gp << "set xlabel 'Estimated diameter D'\n";
    gp << "set ylabel 'N(>=D) / surface area'\n";
    gp << "set logscale xy\n";
    gp << "set grid\n";
    gp << "plot '" << output_dir
       << "/boulder_size_cumulative.csv' using 1:3 with linespoints lw 2 title 'Cumulative density'\n";
}

Eigen::Vector3d ArrayToVec3(const std::array<double, 3>& c) {
    return Eigen::Vector3d(c[0], c[1], c[2]);
}

const std::vector<int>& BoundaryOrAllVertices(const ComponentInfo& info) {
    if (!info.boundary_vertices.empty()) return info.boundary_vertices;
    return info.vertex_indices;
}

double MinDistanceBetweenVertexSets(const std::vector<Eigen::Vector3d>& vertices,
                                    const std::vector<int>& a,
                                    const std::vector<int>& b,
                                    double early_exit_threshold) {
    if (a.empty() || b.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const std::vector<int>* small = &a;
    const std::vector<int>* large = &b;
    if (a.size() > b.size()) std::swap(small, large);

    double best2 = std::numeric_limits<double>::infinity();
    const double early2 = early_exit_threshold * early_exit_threshold;

    for (int ia : *small) {
        const Eigen::Vector3d& pa = vertices[ia];
        for (int ib : *large) {
            const double d2 = (pa - vertices[ib]).squaredNorm();
            if (d2 < best2) {
                best2 = d2;
                if (best2 <= early2) {
                    return std::sqrt(best2);
                }
            }
        }
    }
    return std::sqrt(best2);
}

double MeanFinite(const std::vector<double>& values) {
    double sum = 0.0;
    std::size_t count = 0;
    for (double v : values) {
        if (std::isfinite(v)) {
            sum += v;
            ++count;
        }
    }
    return (count == 0) ? 0.0 : sum / static_cast<double>(count);
}

double MaxFinite(const std::vector<double>& values) {
    double m = -std::numeric_limits<double>::infinity();
    bool found = false;
    for (double v : values) {
        if (std::isfinite(v)) {
            m = std::max(m, v);
            found = true;
        }
    }
    return found ? m : 0.0;
}

double StdFinite(const std::vector<double>& values, double mean) {
    double acc = 0.0;
    std::size_t count = 0;
    for (double v : values) {
        if (std::isfinite(v)) {
            const double d = v - mean;
            acc += d * d;
            ++count;
        }
    }
    return (count == 0) ? 0.0 : std::sqrt(acc / static_cast<double>(count));
}

std::vector<double> GatherFiniteByIndices(const std::vector<double>& src,
                                          const std::vector<int>& indices) {
    std::vector<double> out;
    out.reserve(indices.size());
    for (int idx : indices) {
        if (idx >= 0 && idx < static_cast<int>(src.size()) && std::isfinite(src[idx])) {
            out.push_back(src[idx]);
        }
    }
    return out;
}

std::vector<Eigen::Vector3d> GatherPoints(const std::vector<Eigen::Vector3d>& vertices,
                                          const std::vector<int>& indices) {
    std::vector<Eigen::Vector3d> out;
    out.reserve(indices.size());
    for (int idx : indices) {
        out.push_back(vertices[idx]);
    }
    return out;
}

std::string ToJsonArray(const std::vector<int>& values) {
    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << values[i];
    }
    oss << "]";
    return oss.str();
}

std::string ToJsonArray(const std::vector<double>& values) {
    std::ostringstream oss;
    oss << std::setprecision(16);
    oss << "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << values[i];
    }
    oss << "]";
    return oss.str();
}

std::string ToJsonArray(const std::array<double, 3>& values) {
    std::ostringstream oss;
    oss << std::setprecision(16);
    oss << "[" << values[0] << ", " << values[1] << ", " << values[2] << "]";
    return oss.str();
}

void WriteComponentInfoJson(std::ostream& os, const ComponentInfo& info, int indent) {
    const std::string pad(indent, ' ');
    const std::string pad2(indent + 2, ' ');
    os << pad << "{\n";
    os << pad2 << "\"component_id\": " << info.component_id << ",\n";
    os << pad2 << "\"vertex_indices\": " << ToJsonArray(info.vertex_indices) << ",\n";
    os << pad2 << "\"num_vertices\": " << info.num_vertices << ",\n";
    os << pad2 << "\"centroid\": " << ToJsonArray(info.centroid) << ",\n";
    os << pad2 << "\"max_height\": " << info.max_height << ",\n";
    os << pad2 << "\"mean_height\": " << info.mean_height << ",\n";
    os << pad2 << "\"diameter_estimate\": " << info.diameter_estimate << ",\n";
    os << pad2 << "\"horizontal_size\": " << info.horizontal_size << ",\n";
    os << pad2 << "\"surface_area_estimate\": " << info.surface_area_estimate << ",\n";
    os << pad2 << "\"compactness_ratio\": " << info.compactness_ratio << ",\n";
    os << pad2 << "\"elongation\": " << info.elongation << ",\n";
    os << pad2 << "\"pca_eigenvalues\": " << ToJsonArray(info.pca_eigenvalues) << ",\n";
    os << pad2 << "\"inside_mean_height\": " << info.inside_mean_height << ",\n";
    os << pad2 << "\"inside_max_height\": " << info.inside_max_height << ",\n";
    os << pad2 << "\"boundary_mean_height\": " << info.boundary_mean_height << ",\n";
    os << pad2 << "\"boundary_height_std\": " << info.boundary_height_std << ",\n";
    os << pad2 << "\"outer_ring_mean_height\": " << info.outer_ring_mean_height << ",\n";
    os << pad2 << "\"outer_ring_height_std\": " << info.outer_ring_height_std << ",\n";
    os << pad2 << "\"prominence_mean\": " << info.prominence_mean << ",\n";
    os << pad2 << "\"prominence_max\": " << info.prominence_max << ",\n";
    os << pad2 << "\"boundary_drop\": " << info.boundary_drop << ",\n";
    os << pad2 << "\"boundary_irregularity\": " << info.boundary_irregularity << ",\n";
    os << pad2 << "\"boundary_size\": " << info.boundary_size << ",\n";
    os << pad2 << "\"outer_ring_size\": " << info.outer_ring_size << ",\n";
    os << pad2 << "\"boundary_vertices\": " << ToJsonArray(info.boundary_vertices) << ",\n";
    os << pad2 << "\"outer_ring_vertices\": " << ToJsonArray(info.outer_ring_vertices) << "\n";
    os << pad << "}";
}

std::vector<Eigen::Vector3d> GetVertices(const TriangleMesh& mesh) {
    return mesh.vertices_;
}

std::vector<Eigen::Vector3i> GetTriangles(const TriangleMesh& mesh) {
    return mesh.triangles_;
}

void VisualizeMesh(const std::shared_ptr<TriangleMesh>& mesh, const std::string& title) {
    std::vector<std::shared_ptr<const open3d::geometry::Geometry>> geoms;
    geoms.push_back(mesh);
    open3d::visualization::DrawGeometries(geoms, title);
}

}  // namespace

std::shared_ptr<TriangleMesh> LoadMesh(const std::string& mesh_path, bool verbose) {
    auto mesh = std::make_shared<TriangleMesh>();
    if (!open3d::io::ReadTriangleMesh(mesh_path, *mesh)) {
        throw std::runtime_error("Cannot load mesh from: " + mesh_path);
    }
    if (mesh->IsEmpty()) {
        throw std::runtime_error("Loaded mesh is empty: " + mesh_path);
    }

    mesh->RemoveDuplicatedVertices();
    mesh->RemoveDegenerateTriangles();
    mesh->RemoveDuplicatedTriangles();
    mesh->RemoveUnreferencedVertices();
    mesh->ComputeVertexNormals();
    mesh->ComputeTriangleNormals();

    if (verbose) {
        std::cout << "[INFO] Mesh loaded: " << mesh_path << "\n";
        std::cout << "[INFO] #Vertices = " << mesh->vertices_.size() << "\n";
        std::cout << "[INFO] #Faces    = " << mesh->triangles_.size() << "\n";
    }
    return mesh;
}

std::shared_ptr<TriangleMesh> OptionalSmoothMesh(const TriangleMesh& mesh, int iterations) {
    if (iterations <= 0) {
        return std::make_shared<TriangleMesh>(mesh);
    }
    auto smoothed = mesh.FilterSmoothSimple(iterations);
    smoothed->ComputeVertexNormals();
    smoothed->ComputeTriangleNormals();
    return smoothed;
}

double ComputeMeanEdgeLength(const TriangleMesh& mesh) {
    const auto& vertices = mesh.vertices_;
    const auto& triangles = mesh.triangles_;

    std::set<std::pair<int, int>> edge_set;
    for (const auto& tri : triangles) {
        const int i = tri(0);
        const int j = tri(1);
        const int k = tri(2);
        edge_set.insert({std::min(i, j), std::max(i, j)});
        edge_set.insert({std::min(j, k), std::max(j, k)});
        edge_set.insert({std::min(k, i), std::max(k, i)});
    }

    if (edge_set.empty()) {
        return 0.0;
    }

    double sum = 0.0;
    for (const auto& e : edge_set) {
        sum += (vertices[e.first] - vertices[e.second]).norm();
    }
    return sum / static_cast<double>(edge_set.size());
}

std::vector<std::vector<int>> BuildTopologyAdjacency(const TriangleMesh& mesh) {
    const auto n = static_cast<int>(mesh.vertices_.size());
    std::vector<std::unordered_set<int>> adjacency_set(n);

    for (const auto& tri : mesh.triangles_) {
        const int i = tri(0);
        const int j = tri(1);
        const int k = tri(2);
        adjacency_set[i].insert(j);
        adjacency_set[i].insert(k);
        adjacency_set[j].insert(i);
        adjacency_set[j].insert(k);
        adjacency_set[k].insert(i);
        adjacency_set[k].insert(j);
    }

    std::vector<std::vector<int>> adjacency(n);
    for (int i = 0; i < n; ++i) {
        adjacency[i].assign(adjacency_set[i].begin(), adjacency_set[i].end());
    }
    return adjacency;
}

std::vector<double> TriangleAreas(const TriangleMesh& mesh) {
    std::vector<double> areas(mesh.triangles_.size(), 0.0);
    for (std::size_t t = 0; t < mesh.triangles_.size(); ++t) {
        const auto& tri = mesh.triangles_[t];
        const Eigen::Vector3d& a = mesh.vertices_[tri(0)];
        const Eigen::Vector3d& b = mesh.vertices_[tri(1)];
        const Eigen::Vector3d& c = mesh.vertices_[tri(2)];
        areas[t] = 0.5 * ((b - a).cross(c - a)).norm();
    }
    return areas;
}

std::vector<std::vector<int>> BuildVertexToTriangleMap(const TriangleMesh& mesh) {
    std::vector<std::vector<int>> vertex_to_triangles(mesh.vertices_.size());
    for (std::size_t t = 0; t < mesh.triangles_.size(); ++t) {
        const auto& tri = mesh.triangles_[t];
        vertex_to_triangles[tri(0)].push_back(static_cast<int>(t));
        vertex_to_triangles[tri(1)].push_back(static_cast<int>(t));
        vertex_to_triangles[tri(2)].push_back(static_cast<int>(t));
    }
    return vertex_to_triangles;
}

PlaneFitResult FitPlanePCA(const std::vector<Eigen::Vector3d>& points) {
    PlaneFitResult result;
    if (points.empty()) {
        return result;
    }

    for (const auto& p : points) {
        result.centroid += p;
    }
    result.centroid /= static_cast<double>(points.size());

    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (const auto& p : points) {
        const Eigen::Vector3d d = p - result.centroid;
        cov += d * d.transpose();
    }
    cov /= std::max<std::size_t>(points.size(), 1);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
    if (solver.info() != Eigen::Success) {
        return result;
    }

    result.evals = solver.eigenvalues();      // ascending
    result.evecs = solver.eigenvectors();     // columns correspond to evals
    result.normal = result.evecs.col(0).normalized();
    return result;
}

LocalHeightResult ComputeLocalHeights(const TriangleMesh& mesh,
                                      double radius,
                                      int min_neighbors,
                                      bool use_absolute_height,
                                      bool verbose) {
    const auto n = mesh.vertices_.size();
    LocalHeightResult out;
    out.heights.assign(n, std::numeric_limits<double>::quiet_NaN());
    out.valid_mask.assign(n, 0);
    out.local_roughness.assign(n, std::numeric_limits<double>::quiet_NaN());
    out.local_planarity.assign(n, std::numeric_limits<double>::quiet_NaN());
    out.neighbor_counts.assign(n, 0);

    PointCloud pcd;
    pcd.points_ = mesh.vertices_;
    KDTreeFlann tree;
    tree.SetGeometry(pcd);

    for (std::size_t i = 0; i < n; ++i) {
        const Eigen::Vector3d& p = mesh.vertices_[i];

        std::vector<int> idx;
        std::vector<double> dist2;
        tree.SearchRadius(p, radius, idx, dist2);
        out.neighbor_counts[i] = static_cast<int>(idx.size());

        if (static_cast<int>(idx.size()) < min_neighbors) {
            continue;
        }

        std::vector<Eigen::Vector3d> nbr_points;
        nbr_points.reserve(idx.size());
        for (int id : idx) {
            nbr_points.push_back(mesh.vertices_[id]);
        }

        const PlaneFitResult fit = FitPlanePCA(nbr_points);
        Eigen::Vector3d normal = fit.normal;
        const Eigen::Vector3d& v_normal = mesh.vertex_normals_[i];
        if (normal.dot(v_normal) < 0.0) {
            normal = -normal;
        }

        const double signed_height = normal.dot(p - fit.centroid);
        std::vector<double> signed_distances;
        signed_distances.reserve(nbr_points.size());
        for (const auto& q : nbr_points) {
            signed_distances.push_back((q - fit.centroid).dot(normal));
        }

        out.heights[i] = use_absolute_height ? std::abs(signed_height) : signed_height;
        const double mean_sd = MeanFinite(signed_distances);
        out.local_roughness[i] = StdFinite(signed_distances, mean_sd);
        const double denom = fit.evals.sum() + kEps;
        out.local_planarity[i] = fit.evals(0) / denom;
        out.valid_mask[i] = 1;

        if (verbose && i > 0 && i % 50000 == 0) {
            std::cout << "[INFO] compute_local_heights(radius=" << radius
                      << "): processed " << i << "/" << n << "\n";
        }
    }

    return out;
}

CandidateResult DetectCandidateVertices(const std::vector<double>& heights,
                                        const std::vector<std::uint8_t>& valid_mask,
                                        const std::vector<double>* local_roughness,
                                        double height_threshold,
                                        double height_sigma_factor,
                                        double roughness_max,
                                        bool verbose) {
    CandidateResult out;
    out.candidate_mask.assign(valid_mask.size(), 0);

    std::vector<double> valid_heights;
    valid_heights.reserve(heights.size());
    for (std::size_t i = 0; i < heights.size(); ++i) {
        if (valid_mask[i] && std::isfinite(heights[i])) {
            valid_heights.push_back(heights[i]);
        }
    }

    if (valid_heights.empty()) {
        throw std::runtime_error("No valid heights available.");
    }

    out.mu = MeanFinite(valid_heights);
    out.sigma = StdFinite(valid_heights, out.mu);
    out.used_threshold = std::isfinite(height_threshold)
                                 ? height_threshold
                                 : out.mu + height_sigma_factor * out.sigma;

    for (std::size_t i = 0; i < heights.size(); ++i) {
        if (!valid_mask[i] || !std::isfinite(heights[i])) {
            continue;
        }
        bool ok = heights[i] > out.used_threshold;
        if (ok && local_roughness != nullptr && std::isfinite(roughness_max)) {
            ok = ((*local_roughness)[i] <= roughness_max);
        }
        out.candidate_mask[i] = ok ? 1 : 0;
    }

    if (verbose) {
        const int count = std::accumulate(out.candidate_mask.begin(), out.candidate_mask.end(), 0);
        std::cout << "[INFO] Height mean   = " << out.mu << "\n";
        std::cout << "[INFO] Height std    = " << out.sigma << "\n";
        std::cout << "[INFO] Threshold     = " << out.used_threshold << "\n";
        std::cout << "[INFO] #Candidates   = " << count << "\n";
    }
    return out;
}

void ExtractCandidateComponents(const std::vector<std::uint8_t>& candidate_mask,
                                const std::vector<std::vector<int>>& adjacency,
                                std::vector<int>& labels,
                                std::vector<std::vector<int>>& components,
                                bool verbose) {
    const int n = static_cast<int>(candidate_mask.size());
    labels.assign(n, -1);
    components.clear();

    std::vector<std::uint8_t> visited(n, 0);
    int comp_id = 0;

    for (int v = 0; v < n; ++v) {
        if (!candidate_mask[v] || visited[v]) {
            continue;
        }

        std::queue<int> q;
        q.push(v);
        visited[v] = 1;
        components.push_back({});

        while (!q.empty()) {
            const int cur = q.front();
            q.pop();
            labels[cur] = comp_id;
            components.back().push_back(cur);

            for (int nbr : adjacency[cur]) {
                if (candidate_mask[nbr] && !visited[nbr]) {
                    visited[nbr] = 1;
                    q.push(nbr);
                }
            }
        }
        ++comp_id;
    }

    if (verbose) {
        if (components.empty()) {
            std::cout << "[INFO] No candidate vertices found.\n";
        } else {
            std::cout << "[INFO] #Connected components = " << components.size() << "\n";
        }
    }
}

double EstimateComponentDiameter(const std::vector<Eigen::Vector3d>& points) {
    if (points.empty()) {
        return 0.0;
    }
    Eigen::Vector3d pmin = points.front();
    Eigen::Vector3d pmax = points.front();
    for (const auto& p : points) {
        pmin = pmin.cwiseMin(p);
        pmax = pmax.cwiseMax(p);
    }
    return (pmax - pmin).norm();
}

double EstimateComponentHorizontalSize(const std::vector<Eigen::Vector3d>& points) {
    if (points.size() < 3) {
        return 0.0;
    }

    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    for (const auto& p : points) centroid += p;
    centroid /= static_cast<double>(points.size());

    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (const auto& p : points) {
        Eigen::Vector3d d = p - centroid;
        cov += d * d.transpose();
    }
    cov /= std::max<std::size_t>(points.size(), 1);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
    if (solver.info() != Eigen::Success) {
        return 0.0;
    }

    Eigen::Vector3d evals = solver.eigenvalues();
    Eigen::Matrix3d evecs = solver.eigenvectors();
    std::array<int, 3> order = {2, 1, 0};  // descending for self-adjoint solver

    double min1 = std::numeric_limits<double>::infinity();
    double max1 = -std::numeric_limits<double>::infinity();
    double min2 = std::numeric_limits<double>::infinity();
    double max2 = -std::numeric_limits<double>::infinity();

    for (const auto& p : points) {
        Eigen::Vector3d d = p - centroid;
        double proj1 = d.dot(evecs.col(order[0]));
        double proj2 = d.dot(evecs.col(order[1]));
        min1 = std::min(min1, proj1);
        max1 = std::max(max1, proj1);
        min2 = std::min(min2, proj2);
        max2 = std::max(max2, proj2);
    }

    (void)evals;
    return std::max(max1 - min1, max2 - min2);
}

double ComputeComponentArea(const std::vector<int>& component_vertices,
                            const std::vector<std::vector<int>>& vertex_to_triangles,
                            const std::vector<double>& tri_areas) {
    std::unordered_set<int> tri_ids;
    for (int v : component_vertices) {
        for (int tid : vertex_to_triangles[v]) {
            tri_ids.insert(tid);
        }
    }

    double area = 0.0;
    for (int tid : tri_ids) {
        area += tri_areas[tid];
    }
    return area;
}

std::pair<double, std::vector<double>> ComputeComponentPCAShape(
        const std::vector<Eigen::Vector3d>& points) {
    if (points.size() < 3) {
        return {1.0, {0.0, 0.0, 0.0}};
    }

    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    for (const auto& p : points) centroid += p;
    centroid /= static_cast<double>(points.size());

    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (const auto& p : points) {
        const Eigen::Vector3d d = p - centroid;
        cov += d * d.transpose();
    }
    cov /= std::max<std::size_t>(points.size(), 1);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
    if (solver.info() != Eigen::Success) {
        return {1.0, {0.0, 0.0, 0.0}};
    }

    Eigen::Vector3d evals = solver.eigenvalues();  // ascending
    const double l1 = std::max(evals(2), kEps);
    const double l2 = std::max(evals(1), kEps);
    const double l3 = std::max(evals(0), kEps);
    const double elongation = std::sqrt(l1 / l2);
    return {elongation, {l1, l2, l3}};
}

void GetComponentBoundaryAndOuterRing(const std::vector<int>& component_vertices,
                                      const std::vector<std::vector<int>>& adjacency,
                                      int outer_ring_depth,
                                      std::vector<int>& boundary,
                                      std::vector<int>& outer_ring) {
    std::unordered_set<int> comp_set(component_vertices.begin(), component_vertices.end());
    std::unordered_set<int> boundary_set;
    std::unordered_set<int> outer_ring_set;

    for (int v : component_vertices) {
        bool has_outside_neighbor = false;
        for (int nbr : adjacency[v]) {
            if (comp_set.find(nbr) == comp_set.end()) {
                has_outside_neighbor = true;
                outer_ring_set.insert(nbr);
            }
        }
        if (has_outside_neighbor) {
            boundary_set.insert(v);
        }
    }

    if (outer_ring_depth > 1) {
        std::unordered_set<int> frontier = outer_ring_set;
        std::unordered_set<int> visited = outer_ring_set;
        for (int d = 0; d < outer_ring_depth - 1; ++d) {
            std::unordered_set<int> new_frontier;
            for (int v : frontier) {
                for (int nbr : adjacency[v]) {
                    if (comp_set.find(nbr) == comp_set.end() && visited.find(nbr) == visited.end()) {
                        new_frontier.insert(nbr);
                    }
                }
            }
            visited.insert(new_frontier.begin(), new_frontier.end());
            frontier = std::move(new_frontier);
        }
        outer_ring_set = std::move(visited);
    }

    boundary.assign(boundary_set.begin(), boundary_set.end());
    outer_ring.assign(outer_ring_set.begin(), outer_ring_set.end());
    std::sort(boundary.begin(), boundary.end());
    std::sort(outer_ring.begin(), outer_ring.end());
}

BoundaryMetrics ComputeBoundaryMetrics(const std::vector<int>& component_vertices,
                                       const std::vector<double>& heights,
                                       const std::vector<std::vector<int>>& adjacency,
                                       int outer_ring_depth) {
    BoundaryMetrics m;
    GetComponentBoundaryAndOuterRing(component_vertices, adjacency, outer_ring_depth,
                                     m.boundary_vertices, m.outer_ring_vertices);

    const std::vector<double> inside_h = GatherFiniteByIndices(heights, component_vertices);
    const std::vector<double> boundary_h = GatherFiniteByIndices(heights, m.boundary_vertices);
    const std::vector<double> outer_h = GatherFiniteByIndices(heights, m.outer_ring_vertices);

    m.inside_mean_height = MeanFinite(inside_h);
    m.inside_max_height = MaxFinite(inside_h);
    m.boundary_mean_height = MeanFinite(boundary_h);
    m.boundary_height_std = StdFinite(boundary_h, m.boundary_mean_height);
    m.outer_ring_mean_height = MeanFinite(outer_h);
    m.outer_ring_height_std = StdFinite(outer_h, m.outer_ring_mean_height);

    m.prominence_mean = m.inside_mean_height - m.outer_ring_mean_height;
    m.prominence_max = m.inside_max_height - m.outer_ring_mean_height;
    m.boundary_drop = m.boundary_mean_height - m.outer_ring_mean_height;
    m.boundary_irregularity = m.boundary_height_std;

    m.boundary_size = static_cast<int>(m.boundary_vertices.size());
    m.outer_ring_size = static_cast<int>(m.outer_ring_vertices.size());
    return m;
}

std::vector<ComponentInfo> MeasureComponents(const TriangleMesh& mesh,
                                             const std::vector<std::vector<int>>& components,
                                             const std::vector<double>& heights,
                                             const std::vector<std::vector<int>>& adjacency,
                                             const std::vector<std::vector<int>>& vertex_to_triangles,
                                             const std::vector<double>& tri_areas,
                                             int outer_ring_depth,
                                             bool verbose) {
    std::vector<ComponentInfo> infos;
    infos.reserve(components.size());

    for (std::size_t comp_id = 0; comp_id < components.size(); ++comp_id) {
        const auto& comp_vertices = components[comp_id];
        const std::vector<Eigen::Vector3d> pts = GatherPoints(mesh.vertices_, comp_vertices);
        const std::vector<double> h = GatherFiniteByIndices(heights, comp_vertices);

        if (pts.empty() || h.empty()) {
            continue;
        }

        Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
        for (const auto& p : pts) centroid += p;
        centroid /= static_cast<double>(pts.size());

        ComponentInfo info;
        info.component_id = static_cast<int>(comp_id);
        info.vertex_indices = comp_vertices;
        info.num_vertices = static_cast<int>(comp_vertices.size());
        info.centroid = {centroid(0), centroid(1), centroid(2)};
        info.max_height = MaxFinite(h);
        info.mean_height = MeanFinite(h);
        info.diameter_estimate = EstimateComponentDiameter(pts);
        info.horizontal_size = EstimateComponentHorizontalSize(pts);
        info.surface_area_estimate = ComputeComponentArea(comp_vertices, vertex_to_triangles, tri_areas);
        info.compactness_ratio = info.max_height / (info.horizontal_size + kEps);

        auto [elongation, pca_evals] = ComputeComponentPCAShape(pts);
        info.elongation = elongation;
        info.pca_eigenvalues = std::move(pca_evals);

        const BoundaryMetrics bm = ComputeBoundaryMetrics(comp_vertices, heights, adjacency, outer_ring_depth);
        info.inside_mean_height = bm.inside_mean_height;
        info.inside_max_height = bm.inside_max_height;
        info.boundary_mean_height = bm.boundary_mean_height;
        info.boundary_height_std = bm.boundary_height_std;
        info.outer_ring_mean_height = bm.outer_ring_mean_height;
        info.outer_ring_height_std = bm.outer_ring_height_std;
        info.prominence_mean = bm.prominence_mean;
        info.prominence_max = bm.prominence_max;
        info.boundary_drop = bm.boundary_drop;
        info.boundary_irregularity = bm.boundary_irregularity;
        info.boundary_size = bm.boundary_size;
        info.outer_ring_size = bm.outer_ring_size;
        info.boundary_vertices = bm.boundary_vertices;
        info.outer_ring_vertices = bm.outer_ring_vertices;

        infos.push_back(std::move(info));
    }

    if (verbose) {
        std::cout << "[INFO] Measured " << infos.size() << " components.\n";
    }
    return infos;
}

std::vector<std::vector<int>> BuildComponentMergeGroups(
        const TriangleMesh& mesh,
        const std::vector<ComponentInfo>& component_infos,
        double mean_edge_length,
        int merge_knn,
        double merge_boundary_distance_factor,
        double merge_centroid_distance_factor,
        double merge_height_difference_factor,
        double merge_prominence_difference_factor,
        bool verbose) {
    std::vector<std::vector<int>> groups;

    const int n = static_cast<int>(component_infos.size());
    if (n == 0) return groups;
    if (n == 1) {
        groups.push_back({0});
        return groups;
    }

    PointCloud centroid_cloud;
    centroid_cloud.points_.reserve(component_infos.size());
    for (const auto& info : component_infos) {
        centroid_cloud.points_.push_back(ArrayToVec3(info.centroid));
    }

    KDTreeFlann centroid_tree;
    centroid_tree.SetGeometry(centroid_cloud);

    DisjointSet dsu(n);

    const double boundary_thresh = merge_boundary_distance_factor * mean_edge_length;
    const double height_thresh = merge_height_difference_factor * mean_edge_length;
    const double prominence_thresh = merge_prominence_difference_factor * mean_edge_length;
    const int knn = std::max(2, merge_knn);

    int merged_pair_count = 0;

    for (int i = 0; i < n; ++i) {
        const auto& a = component_infos[i];
        const Eigen::Vector3d ca = ArrayToVec3(a.centroid);

        std::vector<int> idx(knn);
        std::vector<double> dist2(knn);
        const int found = centroid_tree.SearchKNN(ca, std::min(knn, n), idx, dist2);

        for (int t = 0; t < found; ++t) {
            const int j = idx[t];
            if (j <= i) continue;

            const auto& b = component_infos[j];
            const Eigen::Vector3d cb = ArrayToVec3(b.centroid);

            const double centroid_dist = (ca - cb).norm();
            const double local_size =
                    std::max({a.diameter_estimate, b.diameter_estimate, 3.0 * mean_edge_length});
            const double centroid_thresh = merge_centroid_distance_factor * local_size;

            if (centroid_dist > centroid_thresh) continue;
            if (std::abs(a.mean_height - b.mean_height) > height_thresh) continue;
            if (std::abs(a.prominence_mean - b.prominence_mean) > prominence_thresh) continue;

            const auto& va = BoundaryOrAllVertices(a);
            const auto& vb = BoundaryOrAllVertices(b);

            const double boundary_dist = MinDistanceBetweenVertexSets(
                    mesh.vertices_, va, vb, boundary_thresh);

            if (boundary_dist > boundary_thresh) continue;

            dsu.Union(i, j);
            ++merged_pair_count;
        }
    }

    std::unordered_map<int, std::vector<int>> root_to_group;
    for (int i = 0; i < n; ++i) {
        root_to_group[dsu.Find(i)].push_back(i);
    }

    groups.reserve(root_to_group.size());
    for (auto& kv : root_to_group) {
        groups.push_back(std::move(kv.second));
    }

    std::sort(groups.begin(), groups.end(),
              [](const std::vector<int>& a, const std::vector<int>& b) {
                  return a.front() < b.front();
              });

    if (verbose) {
        int multi_groups = 0;
        for (const auto& g : groups) {
            if (g.size() > 1) ++multi_groups;
        }
        std::cout << "[INFO] Component merge pairs   = " << merged_pair_count << "\n";
        std::cout << "[INFO] Merge groups total      = " << groups.size() << "\n";
        std::cout << "[INFO] Merge groups (size > 1) = " << multi_groups << "\n";
    }

    return groups;
}

std::vector<std::vector<int>> BuildMergedComponents(
        const std::vector<ComponentInfo>& component_infos,
        const std::vector<std::vector<int>>& merge_groups) {
    std::vector<std::vector<int>> merged_components;
    merged_components.reserve(merge_groups.size());

    for (const auto& group : merge_groups) {
        std::unordered_set<int> vertex_set;

        for (int comp_idx : group) {
            const auto& verts = component_infos[comp_idx].vertex_indices;
            vertex_set.insert(verts.begin(), verts.end());
        }

        std::vector<int> merged(vertex_set.begin(), vertex_set.end());
        std::sort(merged.begin(), merged.end());
        merged_components.push_back(std::move(merged));
    }

    return merged_components;
}
void FilterBoulderComponents(const std::vector<ComponentInfo>& component_infos,
                             std::vector<ComponentInfo>& kept,
                             std::vector<ComponentInfo>& rejected,
                             int min_vertices,
                             double min_max_height,
                             double min_diameter,
                             double max_diameter,
                             double min_compactness,
                             double max_elongation,
                             double min_prominence_mean,
                             double min_prominence_max,
                             double min_boundary_drop,
                             double min_boundary_irregularity,
                             double max_boundary_irregularity,
                             bool verbose) {
    kept.clear();
    rejected.clear();

    for (const auto& info : component_infos) {
        bool ok = true;
        if (info.num_vertices < min_vertices) ok = false;
        if (info.max_height < min_max_height) ok = false;
        if (info.diameter_estimate < min_diameter) ok = false;
        if (info.diameter_estimate > max_diameter) ok = false;
        if (info.compactness_ratio < min_compactness) ok = false;
        if (info.elongation > max_elongation) ok = false;
        if (info.prominence_mean < min_prominence_mean) ok = false;
        if (info.prominence_max < min_prominence_max) ok = false;
        if (info.boundary_drop < min_boundary_drop) ok = false;
        if (info.boundary_irregularity < min_boundary_irregularity) ok = false;
        if (info.boundary_irregularity > max_boundary_irregularity) ok = false;

        if (ok) kept.push_back(info);
        else rejected.push_back(info);
    }

    if (verbose) {
        std::cout << "[INFO] Kept components     = " << kept.size() << "\n";
        std::cout << "[INFO] Rejected components = " << rejected.size() << "\n";
    }
}

BoulderSizeDistribution ComputeBoulderSizeDistribution(
        const std::vector<ComponentInfo>& kept_boulders,
        double total_surface_area,
        int num_bins,
        bool log_bins,
        bool verbose) {
    BoulderSizeDistribution out;

    for (const auto& info : kept_boulders) {
        if (std::isfinite(info.diameter_estimate) && info.diameter_estimate > 0.0) {
            out.diameters_sorted.push_back(info.diameter_estimate);
        }
    }

    if (out.diameters_sorted.empty()) {
        if (verbose) {
            std::cout << "[INFO] No valid boulder diameters for size distribution.\n";
        }
        return out;
    }

    std::sort(out.diameters_sorted.begin(), out.diameters_sorted.end());

    const int bins = std::max(1, num_bins);
    const double dmin = out.diameters_sorted.front();
    const double dmax = out.diameters_sorted.back();

    out.bin_edges.resize(bins + 1, dmin);
    out.bin_centers.resize(bins, dmin);
    out.bin_counts.assign(bins, 0);
    out.cumulative_counts_ge.assign(bins, 0);
    out.cumulative_number_density_ge.assign(bins, 0.0);

    if (dmax <= dmin * (1.0 + 1e-12)) {
        out.bin_edges[0] = dmin;
        out.bin_edges[1] = dmax * (1.0 + 1e-6) + 1e-12;
        out.bin_centers[0] = dmin;
        out.bin_counts[0] = static_cast<int>(out.diameters_sorted.size());
        out.cumulative_counts_ge[0] = static_cast<int>(out.diameters_sorted.size());
        out.cumulative_number_density_ge[0] =
                static_cast<double>(out.cumulative_counts_ge[0]) /
                std::max(total_surface_area, 1e-12);
        return out;
    }

    if (log_bins && dmin > 0.0) {
        const double log_min = std::log10(dmin);
        const double log_max = std::log10(dmax);

        for (int i = 0; i <= bins; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(bins);
            out.bin_edges[i] = std::pow(10.0, log_min + t * (log_max - log_min));
        }
        for (int i = 0; i < bins; ++i) {
            out.bin_centers[i] = std::sqrt(out.bin_edges[i] * out.bin_edges[i + 1]);
        }
    } else {
        const double step = (dmax - dmin) / static_cast<double>(bins);
        for (int i = 0; i <= bins; ++i) {
            out.bin_edges[i] = dmin + step * static_cast<double>(i);
        }
        for (int i = 0; i < bins; ++i) {
            out.bin_centers[i] = 0.5 * (out.bin_edges[i] + out.bin_edges[i + 1]);
        }
    }

    for (double d : out.diameters_sorted) {
        auto it = std::upper_bound(out.bin_edges.begin(), out.bin_edges.end(), d);
        int idx = static_cast<int>(it - out.bin_edges.begin()) - 1;
        idx = std::max(0, std::min(idx, bins - 1));
        out.bin_counts[idx] += 1;
    }

    int running = 0;
    for (int i = bins - 1; i >= 0; --i) {
        running += out.bin_counts[i];
        out.cumulative_counts_ge[i] = running;
        out.cumulative_number_density_ge[i] =
                static_cast<double>(running) / std::max(total_surface_area, 1e-12);
    }

    if (verbose) {
        std::cout << "[INFO] Boulder size distribution computed for "
                  << out.diameters_sorted.size() << " boulders.\n";
    }

    return out;
}

std::vector<Eigen::Vector3d> ScalarToColors(const std::vector<double>& values,
                                            const std::vector<std::uint8_t>* valid_mask) {
    std::vector<Eigen::Vector3d> colors(values.size(), Eigen::Vector3d(0.0, 0.0, 0.0));

    std::vector<std::uint8_t> local_mask;
    if (valid_mask == nullptr) {
        local_mask.resize(values.size(), 0);
        for (std::size_t i = 0; i < values.size(); ++i) {
            local_mask[i] = std::isfinite(values[i]) ? 1 : 0;
        }
        valid_mask = &local_mask;
    }

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    bool found = false;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if ((*valid_mask)[i] && std::isfinite(values[i])) {
            vmin = std::min(vmin, values[i]);
            vmax = std::max(vmax, values[i]);
            found = true;
        }
    }

    if (!found) return colors;

    if (std::abs(vmax - vmin) < kEps) {
        for (std::size_t i = 0; i < values.size(); ++i) {
            if ((*valid_mask)[i]) colors[i] = Eigen::Vector3d(0.5, 0.5, 0.5);
        }
        return colors;
    }

    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!(*valid_mask)[i] || !std::isfinite(values[i])) {
            colors[i] = Eigen::Vector3d(0.2, 0.2, 0.2);
            continue;
        }

        double x = (values[i] - vmin) / (vmax - vmin);
        x = std::clamp(x, 0.0, 1.0);
        const double r = std::clamp(1.5 - std::abs(4.0 * x - 3.0), 0.0, 1.0);
        const double g = std::clamp(1.5 - std::abs(4.0 * x - 2.0), 0.0, 1.0);
        const double b = std::clamp(1.5 - std::abs(4.0 * x - 1.0), 0.0, 1.0);
        colors[i] = Eigen::Vector3d(r, g, b);
    }
    return colors;
}

std::shared_ptr<TriangleMesh> ColorizeScalarField(const TriangleMesh& mesh,
                                                  const std::vector<double>& scalars,
                                                  const std::vector<std::uint8_t>* valid_mask) {
    auto out = std::make_shared<TriangleMesh>(mesh);
    out->vertex_colors_ = ScalarToColors(scalars, valid_mask);
    return out;
}

std::shared_ptr<TriangleMesh> ColorizeCandidateVertices(const TriangleMesh& mesh,
                                                        const std::vector<std::uint8_t>& candidate_mask,
                                                        const Eigen::Vector3d& base_color,
                                                        const Eigen::Vector3d& candidate_color) {
    auto out = std::make_shared<TriangleMesh>(mesh);
    out->vertex_colors_.assign(mesh.vertices_.size(), base_color);
    for (std::size_t i = 0; i < candidate_mask.size(); ++i) {
        if (candidate_mask[i]) {
            out->vertex_colors_[i] = candidate_color;
        }
    }
    return out;
}

std::shared_ptr<TriangleMesh> ColorizeComponents(const TriangleMesh& mesh,
                                                 const std::vector<ComponentInfo>& component_infos) {
    auto out = std::make_shared<TriangleMesh>(mesh);
    out->vertex_colors_.assign(mesh.vertices_.size(), Eigen::Vector3d(0.65, 0.65, 0.65));

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    for (const auto& info : component_infos) {
        const Eigen::Vector3d color(unif(rng), unif(rng), unif(rng));
        for (int vid : info.vertex_indices) {
            out->vertex_colors_[vid] = color;
        }
    }
    return out;
}

bool SaveMesh(const TriangleMesh& mesh, const std::string& output_path) {
    const bool ok = open3d::io::WriteTriangleMesh(output_path, mesh,
                                                  false, false,
                                                  true, true, true, false);
    if (ok) {
        std::cout << "[INFO] Saved mesh: " << output_path << "\n";
    } else {
        std::cout << "[WARN] Failed to save mesh: " << output_path << "\n";
    }
    return ok;
}

void SaveComponentInfosJson(const std::vector<ComponentInfo>& component_infos,
                            const std::string& output_json_path) {
    std::ofstream ofs(output_json_path);
    ofs << "[\n";
    for (std::size_t i = 0; i < component_infos.size(); ++i) {
        WriteComponentInfoJson(ofs, component_infos[i], 2);
        if (i + 1 < component_infos.size()) {
            ofs << ",";
        }
        ofs << "\n";
    }
    ofs << "]\n";
}

void SaveDoubleVectorTxt(const std::vector<double>& values, const std::string& path) {
    std::ofstream ofs(path);
    ofs << std::setprecision(16);
    for (double v : values) {
        ofs << v << '\n';
    }
}

void SaveIntVectorTxt(const std::vector<int>& values, const std::string& path) {
    std::ofstream ofs(path);
    for (int v : values) {
        ofs << v << '\n';
    }
}

void SaveMaskTxt(const std::vector<std::uint8_t>& values, const std::string& path) {
    std::ofstream ofs(path);
    for (std::uint8_t v : values) {
        ofs << static_cast<int>(v) << '\n';
    }
}

DetectionResult DetectBouldersSingleScale(const PipelineParams& params) {
    if (params.mesh_path.empty()) {
        throw std::runtime_error("PipelineParams.mesh_path must not be empty.");
    }

    std::filesystem::create_directories(params.output_dir);

    auto mesh = LoadMesh(params.mesh_path, params.verbose);
    mesh = OptionalSmoothMesh(*mesh, params.smooth_iterations);
    SaveMesh(*mesh, params.output_dir + "/mesh_cleaned.ply");

    const double mean_edge = ComputeMeanEdgeLength(*mesh);
    const double radius = params.radius_multiplier * mean_edge;

    std::cout << "[INFO] Mean edge length = " << mean_edge << "\n";
    std::cout << "[INFO] Radius = " << radius << " ("
              << params.radius_multiplier << " * mean_edge_length)\n";

    const auto adjacency = BuildTopologyAdjacency(*mesh);
    const auto tri_areas = TriangleAreas(*mesh);
    const auto vertex_to_triangles = BuildVertexToTriangleMap(*mesh);

    const LocalHeightResult local = ComputeLocalHeights(*mesh, radius,
                                                        params.min_neighbors,
                                                        params.use_absolute_height,
                                                        params.verbose);

    SaveDoubleVectorTxt(local.heights, params.output_dir + "/heights.txt");
    SaveMaskTxt(local.valid_mask, params.output_dir + "/valid_mask.txt");
    SaveDoubleVectorTxt(local.local_roughness, params.output_dir + "/local_roughness.txt");
    SaveDoubleVectorTxt(local.local_planarity, params.output_dir + "/local_planarity.txt");
    SaveIntVectorTxt(local.neighbor_counts, params.output_dir + "/neighbor_counts.txt");

    const CandidateResult cand = DetectCandidateVertices(local.heights,
                                                         local.valid_mask,
                                                         &local.local_roughness,
                                                         params.height_threshold,
                                                         params.height_sigma_factor,
                                                         params.roughness_max,
                                                         params.verbose);
    SaveMaskTxt(cand.candidate_mask, params.output_dir + "/candidate_mask.txt");

    std::vector<int> labels;
    std::vector<std::vector<int>> components;
    ExtractCandidateComponents(cand.candidate_mask, adjacency, labels, components, params.verbose);
    SaveIntVectorTxt(labels, params.output_dir + "/component_labels.txt");

    
    const auto raw_component_infos = MeasureComponents(*mesh, components, local.heights,
                                                   adjacency, vertex_to_triangles,
                                                   tri_areas, params.outer_ring_depth,
                                                   params.verbose);

    SaveComponentInfosJson(raw_component_infos, params.output_dir + "/all_components_raw.json");

    std::vector<std::vector<int>> working_components = components;
    std::vector<ComponentInfo> working_component_infos = raw_component_infos;

    if (params.enable_component_merging) {
        const auto merge_groups = BuildComponentMergeGroups(
                *mesh,
                raw_component_infos,
                mean_edge,
                params.merge_knn,
                params.merge_boundary_distance_factor,
                params.merge_centroid_distance_factor,
                params.merge_height_difference_factor,
                params.merge_prominence_difference_factor,
                params.verbose);

        working_components = BuildMergedComponents(raw_component_infos, merge_groups);

        working_component_infos = MeasureComponents(*mesh, working_components, local.heights,
                                                adjacency, vertex_to_triangles,
                                                tri_areas, params.outer_ring_depth,
                                                params.verbose);

        SaveComponentInfosJson(working_component_infos,
                            params.output_dir + "/all_components_merged.json");
        SaveIntVectorTxt(std::vector<int>{static_cast<int>(working_components.size())},
                        params.output_dir + "/merged_component_count.txt");
    } else {
        SaveComponentInfosJson(working_component_infos,
                            params.output_dir + "/all_components_merged.json");
    }

    std::vector<ComponentInfo> kept;
    std::vector<ComponentInfo> rejected;
    FilterBoulderComponents(working_component_infos,
                            kept,
                            rejected,
                            params.min_vertices,
                            params.min_max_height,
                            params.min_diameter,
                            params.max_diameter,
                            params.min_compactness,
                            params.max_elongation,
                            params.min_prominence_mean,
                            params.min_prominence_max,
                            params.min_boundary_drop,
                            params.min_boundary_irregularity,
                            params.max_boundary_irregularity,
                            params.verbose);

    SaveComponentInfosJson(kept, params.output_dir + "/kept_boulders.json");
    SaveComponentInfosJson(rejected, params.output_dir + "/rejected_components.json");

    std::vector<std::uint8_t> pred_vertex_mask(mesh->vertices_.size(), 0);
    for (const auto& info : kept) {
        for (int vid : info.vertex_indices) {
            pred_vertex_mask[vid] = 1;
        }
    }
    SaveMaskTxt(pred_vertex_mask, params.output_dir + "/pred_vertex_mask.txt");

    auto height_mesh = ColorizeScalarField(*mesh, local.heights, &local.valid_mask);
    auto candidate_mesh = ColorizeCandidateVertices(*mesh, cand.candidate_mask);
    auto kept_mesh = ColorizeComponents(*mesh, kept);

    SaveMesh(*height_mesh, params.output_dir + "/height_field.ply");
    SaveMesh(*candidate_mesh, params.output_dir + "/candidate_vertices.ply");
    SaveMesh(*kept_mesh, params.output_dir + "/detected_boulders.ply");
    const double total_surface_area =
        std::accumulate(tri_areas.begin(), tri_areas.end(), 0.0);

    const auto size_dist = ComputeBoulderSizeDistribution(
            kept,
            total_surface_area,
            20,     // số bin
            true,   // log bins
            params.verbose);

    SaveDoubleVectorTxt(size_dist.diameters_sorted,
                        params.output_dir + "/boulder_diameters.txt");

    SaveCsvXY_Int(size_dist.bin_centers,
                size_dist.bin_counts,
                params.output_dir + "/boulder_size_histogram.csv",
                "diameter_center",
                "count");

    SaveCsvXY_Int(size_dist.bin_centers,
                size_dist.cumulative_counts_ge,
                params.output_dir + "/boulder_size_cumulative_counts.csv",
                "diameter_center",
                "N_ge_D");

    SaveCsvXY_Double(size_dist.bin_centers,
                    size_dist.cumulative_number_density_ge,
                    params.output_dir + "/boulder_size_cumulative.csv",
                    "diameter_center",
                    "N_ge_D_per_surface_area");

    WriteBoulderDistributionGnuplotScript(params.output_dir, params.mesh_path);
    const std::string gp_cmd =
        "gnuplot \"" + params.output_dir + "/plot_boulder_distribution.gp\"";
    std::system(gp_cmd.c_str());

    if (params.visualize) {
        std::cout << "[INFO] Visualizing height field...\n";
        VisualizeMesh(height_mesh, "Height Field");

        std::cout << "[INFO] Visualizing candidate vertices...\n";
        VisualizeMesh(candidate_mesh, "Candidate Vertices");

        std::cout << "[INFO] Visualizing filtered boulders...\n";
        VisualizeMesh(kept_mesh, "Filtered Boulder-like Protrusions");
    }

    std::cout << "[INFO] Pipeline done.\n";
    std::cout << "[INFO] Output directory: " << params.output_dir << "\n";

    DetectionResult result;
    result.mesh = mesh;
    result.mean_edge_length = mean_edge;
    result.radius = radius;
    result.heights = local.heights;
    result.valid_mask = local.valid_mask;
    result.candidate_mask = cand.candidate_mask;
    result.pred_vertex_mask = pred_vertex_mask;
    result.labels = labels;
    result.components = working_components;
    result.all_component_infos = working_component_infos;
    result.kept_boulders = kept;
    result.rejected_components = rejected;
    result.used_threshold = cand.used_threshold;
    return result;
}

}  // namespace boulder
