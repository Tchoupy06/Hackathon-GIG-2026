#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
//#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>


namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " fichier.obj" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    Mesh mesh;

    if (!CGAL::IO::read_polygon_mesh(filename, mesh)) {
        std::cerr << "Erreur de chargement !" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Mesh chargé : "
        << num_vertices(mesh) << " sommets, "
        << num_faces(mesh) << " faces" << std::endl;

    //PMP::triangulate_faces(mesh);

    double vol = PMP::volume(mesh);
    std::cout << "Volume = " << vol << std::endl;

    double area = PMP::area(mesh);
    std::cout << "Aire = " << area << std::endl;
}