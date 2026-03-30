#pragma once

#include "MeshData.hpp"
#include <iostream>


void loadMesh(std::string filename, Mesh& mesh);

void displayMeshInfo(Mesh& mesh);

void displayMeshData(MeshData& mesh_data);