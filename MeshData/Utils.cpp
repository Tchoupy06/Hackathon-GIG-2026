#include "Utils.hpp"


void loadMesh(std::string filename, Mesh& mesh) {
    if (!CGAL::IO::read_polygon_mesh(filename, mesh)) std::cerr << "Erreur de chargement !" << std::endl;
}

void displayMeshInfo(Mesh& mesh) {
    std::cout << num_vertices(mesh) << " sommets, "
        << num_faces(mesh) << " faces." << std::endl;
}

void displayMeshData(MeshData& mesh_data) {
    std::cout << "Volume = " << mesh_data.volume << " km²"
        << "\nAire = " << mesh_data.area << " km²" << std::endl;
}