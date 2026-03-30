#pragma once

#include <iostream>
#include <fstream>
#include <cstdlib>


struct MeshData {
	double volume = 0;
	double area = 0;
};


MeshData extractMeshData(const std::string& temp_file);


MeshData getMeshData(const std::string& mesh);
