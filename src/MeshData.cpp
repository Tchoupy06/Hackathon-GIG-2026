#include "../include/MeshData.hpp"

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/helpers.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

const char* linearUnitToString(LinearUnit unit) {
    switch (unit) {
        case LinearUnit::Meter:     return "m";
        case LinearUnit::Kilometer: return "km";
        default:                    return "unknown";
    }
}

const char* areaUnitToString(LinearUnit unit) {
    switch (unit) {
        case LinearUnit::Meter:     return "m^2";
        case LinearUnit::Kilometer: return "km^2";
        default:                    return "unknown^2";
    }
}

const char* volumeUnitToString(LinearUnit unit) {
    switch (unit) {
        case LinearUnit::Meter:     return "m^3";
        case LinearUnit::Kilometer: return "km^3";
        default:                    return "unknown^3";
    }
}

double convertVolumeToM3(double volume, LinearUnit unit) {
    switch (unit) {
        case LinearUnit::Meter:
            return volume;
        case LinearUnit::Kilometer:
            return volume * 1e9; // 1 km^3 = 10^9 m^3
        default:
            return volume;
    }
}

double convertAreaToM2(double area, LinearUnit unit) {
    switch (unit) {
        case LinearUnit::Meter:
            return area;
        case LinearUnit::Kilometer:
            return area * 1e6; // 1 km^2 = 10^6 m^2
        default:
            return area;
    }
}

bool loadMesh(const std::string& path, Mesh& mesh) {
    mesh.clear();

    if (!CGAL::IO::read_polygon_mesh(path, mesh)) {
        std::cerr << "Failed to read mesh: " << path << '\n';
        return false;
    }

    if (!CGAL::is_triangle_mesh(mesh)) {
        std::cerr << "Mesh is not triangulated: " << path << '\n';
        return false;
    }

    return true;
}

MeshData getMeshData(const Mesh& mesh, double density, LinearUnit unit) {
    MeshData data;
    data.unit = unit;
    data.density = density;

    data.num_vertices = mesh.number_of_vertices();
    data.num_faces    = mesh.number_of_faces();
    data.num_edges    = mesh.number_of_edges();

    data.is_closed = CGAL::is_closed(mesh);

    // area is always meaningful for a triangle mesh
    data.area = PMP::area(mesh);

    // volume only makes sense if the mesh bounds a volume
    if (data.is_closed && PMP::does_bound_a_volume(mesh)) {
        const double raw_volume = PMP::volume(mesh);

        // In case orientation is inward, keep physical volume positive
        data.volume = std::abs(raw_volume);

        // density is kg/m^3, so convert native volume to m^3 before mass computation
        data.mass = data.density * convertVolumeToM3(data.volume, unit);
        data.has_valid_volume = true;
    } else {
        data.volume = 0.0;
        data.mass = 0.0;
        data.has_valid_volume = false;
    }

    return data;
}

void displayMeshData(const std::string& astro_name, const MeshData& data) {
    std::cout << "====================================\n";
    std::cout << "Astre     : " << astro_name << '\n';
    std::cout << "Unit      : " << linearUnitToString(data.unit) << '\n';
    std::cout << "Vertices  : " << data.num_vertices << '\n';
    std::cout << "Edges     : " << data.num_edges << '\n';
    std::cout << "Faces     : " << data.num_faces << '\n';
    std::cout << "Closed    : " << (data.is_closed ? "yes" : "no") << '\n';

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Area      : " << data.area << ' ' << areaUnitToString(data.unit) << '\n';

    if (data.has_valid_volume) {
        std::cout << "Volume    : " << data.volume << ' ' << volumeUnitToString(data.unit) << '\n';
        std::cout << "Density   : " << data.density << " kg/m^3\n";
        std::cout << "Mass      : " << data.mass << " kg\n";
    } else {
        std::cout << "Volume    : N/A (mesh does not bound a closed volume)\n";
        std::cout << "Density   : " << data.density << " kg/m^3\n";
        std::cout << "Mass      : N/A\n";
    }
}