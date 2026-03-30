#include "MeshData.hpp"


MeshData getMeshData(Mesh& mesh) {
	return MeshData(PMP::volume(mesh), PMP::area(mesh));
}
