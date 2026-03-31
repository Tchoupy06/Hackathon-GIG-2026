#pragma once

#include "MeshTypes.hpp"
#include <string>

struct MeshData {
    double area = 0.0;       // surface area
    double volume = 0.0;     // enclosed volume
    double density = 0.0;    // masse volumique
    double mass = 0.0;       // density * volume

    bool is_closed = false;
    std::size_t num_vertices = 0;
    std::size_t num_faces = 0;
    std::size_t num_edges = 0;
};

bool loadMesh(const std::string& path, Mesh& mesh);
MeshData getMeshData(const Mesh& mesh, double density);
void displayMeshData(const std::string& astro_name, const MeshData& data);