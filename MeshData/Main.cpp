#include "MeshData.hpp"
#include "Utils.hpp"


int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fichier>" << std::endl;
        return 1;
    }

    Mesh mesh;
    loadMesh(argv[1], mesh);

    std::cout << "Maillage chargé : ";
    displayMeshInfo(mesh);

    MeshData mesh_data = getMeshData(mesh);

    displayMeshData(mesh_data);
}