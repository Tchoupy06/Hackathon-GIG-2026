#include "MeshData.hpp"


MeshData extractMeshData(const std::string& temp_file) {
    MeshData mesh_data = MeshData();

    std::ifstream infile(temp_file);

    if (!infile) {
        std::cerr << "Impossible d'ouvrir le fichier " << temp_file << "\n";
        return mesh_data;
    }

    double area, volume;
    infile >> area >> volume; // lit les deux doubles (ignore automatiquement les sauts de ligne)

    if (!infile) {
        std::cerr << "Erreur de lecture des doubles\n";
        return mesh_data;
    }

    infile.close();

    mesh_data.area = area;
    mesh_data.volume = volume;

    return mesh_data;
}


MeshData getMeshData(const std::string& mesh) {
    std::string temp_file = "temp.msd";
    std::string command = "python meshData.py " + mesh + " " + temp_file;

    int ret = std::system(command.c_str()); // ou python3
    if (ret != 0) {
        std::cerr << "Erreur d'exécution Python\n";
        return MeshData();
    }
    
    return extractMeshData(temp_file);
}