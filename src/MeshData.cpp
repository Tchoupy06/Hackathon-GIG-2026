#include "../include/MeshData.hpp"

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/helpers.h>

#include <iostream>
#include <stdexcept>

namespace PMP = CGAL::Polygon_mesh_processing;

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

MeshData getMeshData(const Mesh& mesh, double density) {
    MeshData data;

    data.num_vertices = mesh.number_of_vertices();
    data.num_faces    = mesh.number_of_faces();
    data.num_edges    = mesh.number_of_edges();

    data.is_closed = CGAL::is_closed(mesh);

    data.area = PMP::area(mesh);
    data.density = density;

    if (data.is_closed) {
        data.volume = PMP::volume(mesh);
        data.mass   = data.density * data.volume;
    } else {
        data.volume = 0.0;
        data.mass   = 0.0;
    }

    return data;
}

void displayMeshData(const std::string& astro_name, const MeshData& data) {
    std::cout << "====================================\n";
    std::cout << "Astre     : " << astro_name << '\n';
    std::cout << "Vertices  : " << data.num_vertices << '\n';
    std::cout << "Edges     : " << data.num_edges << '\n';
    std::cout << "Faces     : " << data.num_faces << '\n';
    std::cout << "Closed    : " << (data.is_closed ? "yes" : "no") << '\n';
    std::cout << "Area      : " << data.area << " m^2\n";
    std::cout << "Volume    : " << data.volume << " m^3\n";
    std::cout << "Density   : " << data.density << " kg/m^3\n";
    std::cout << "Mass      : " << data.mass << " kg\n";
}