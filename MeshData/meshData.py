from open3d.cpu.pybind.geometry import TriangleMesh
from open3d.cpu.pybind.io import read_triangle_mesh
from multiprocessing import Pool
from os import cpu_count
from pathlib import Path
from sys import argv
import numpy as np


class Mesh:
    def __init__(self, vertices: np.ndarray, triangles: np.ndarray, surface: float):
        self.vertices: np.ndarray = vertices
        self.triangles: np.ndarray = triangles
        self.surface: float = surface
        self.volume: float = 0

    def chunk_volume(self, args: tuple[int, int]) -> float:
        start_idx, end_idx = args
        volume: float = 0.0
        for i in range(start_idx, end_idx):
            triangle: np.ndarray = self.triangles[i]
            a, b, c = self.vertices[triangle[0]], self.vertices[triangle[1]], self.vertices[triangle[2]]
            volume += np.dot(np.cross(a, b), c) / 6.0

        return volume

    def compute_volume(self, n_processes: int = 4) -> None:
        triangles_number: int = len(self.triangles)
        chunk_size: int = triangles_number // n_processes
        args_list: list[tuple[int, int]] = []

        for i in range(n_processes):
            start: int = i * chunk_size
            end: int = (i + 1) * chunk_size if i != n_processes - 1 else triangles_number
            args_list.append((start, end))

        with Pool(processes=n_processes) as pool:
            results = pool.map(self.chunk_volume, args_list)

        self.volume = abs(sum(results))

    def save(self, filepath: Path) -> None:
        with open(filepath, 'w') as f:
            f.write(f"{self.surface}\n{self.volume}")

    @staticmethod
    def load(filepath: Path) -> 'Mesh':
        triangle_mesh: TriangleMesh = read_triangle_mesh(filepath)
        mesh_vertices: np.ndarray = np.asarray(triangle_mesh.vertices)
        mesh_triangles: np.ndarray = np.asarray(triangle_mesh.triangles)

        return Mesh(mesh_vertices, mesh_triangles, triangle_mesh.get_surface_area())


if __name__ == '__main__':
    if len(argv) != 3:
        print("Usage: python meshData.py <input> <output>")
        exit(1)

    input_mesh_path: Path = Path(argv[1])
    output_mesh_path: Path = Path(argv[2])

    processes_number: int = max(1, cpu_count() - 1)

    mesh: Mesh = Mesh.load(input_mesh_path)
    mesh.compute_volume(processes_number)

    mesh.save(output_mesh_path)
