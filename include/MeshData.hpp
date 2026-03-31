#pragma once

#include "MeshTypes.hpp"

#include <cstddef>
#include <string>

enum class LinearUnit {
    Meter,
    Kilometer
};

struct MeshData {
    double area = 0.0;       // area in native unit^2
    double volume = 0.0;     // volume in native unit^3
    double density = 0.0;    // kg/m^3
    double mass = 0.0;       // kg

    LinearUnit unit = LinearUnit::Meter;

    bool is_closed = false;
    bool has_valid_volume = false;

    std::size_t num_vertices = 0;
    std::size_t num_faces = 0;
    std::size_t num_edges = 0;
};

bool loadMesh(const std::string& path, Mesh& mesh);

MeshData getMeshData(const Mesh& mesh, double density, LinearUnit unit);

void displayMeshData(const std::string& astro_name, const MeshData& data);

const char* linearUnitToString(LinearUnit unit);
const char* areaUnitToString(LinearUnit unit);
const char* volumeUnitToString(LinearUnit unit);

double convertVolumeToM3(double volume, LinearUnit unit);
double convertAreaToM2(double area, LinearUnit unit);