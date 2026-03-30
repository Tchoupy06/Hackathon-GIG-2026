#pragma once

#include "MeshData.hpp"


int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cout << "Usage: MeshData.exe <mesh>" << std::endl;
        return EXIT_FAILURE;
    }

    MeshData mesh_data = getMeshData(argv[1]);

    std::cout << "Aire: " << mesh_data.area << " km2"
        << "\nVolume: " << mesh_data.volume << " km2" << std::endl;

    return EXIT_SUCCESS;
}