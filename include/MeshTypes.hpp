#pragma once

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point  = Kernel::Point_3;
using Mesh   = CGAL::Surface_mesh<Point>;