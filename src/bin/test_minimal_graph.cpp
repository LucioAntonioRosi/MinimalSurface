#include <iostream>
#include <string.h>
#include <functional>

#include "minimal_graph.h"
#include "my_mesh.h"
#include "tiny_blas.h"
#include "P1.h"

static double test_f(const Vec2d &pos)
{
    double x = pos.x;
    double y = pos.y;
    //return 3*x + pow(2*y, 2);
    return 3*x + 900*pow(y, 2);
    //return 5*exp(x) + 1/(1 + y*y) + 2*x*y + pow(5*x,2);
    //return sin(6*M_PI*x) + sin(2*M_PI*y);
    //return 1/(pow(x,2) + pow(y,2));
    //return (pow(x,2) + pow(y,2) >= 0.99 && x >= 0);
    //return 1.0;
}

int main(int argc, char **argv){

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <mesh type> <number of vertices in 1 direction for the square or number of subcircles for the disk>" << std::endl;
        std::cerr << "Mesh type: disk or square" << std::endl;
        return 1;
    }

    const char* mesh_type = argv[1];
    int N = atoi(argv[2]);

    if (N <= 0) {
        std::cerr << "Number of vertices must be a positive integer." << std::endl;
        return 1;
    }

    Mesh m;
    if (strcmp(mesh_type, "disk") == 0) {
        build_disk_mesh(&m, N, 1.0);
    } else if (strcmp(mesh_type, "square") == 0) {
        build_square_mesh(&m, N, 1.0, 1.0);
    } else {
        std::cerr << "Invalid mesh type. Use 'disk' or 'square'." << std::endl;
        return 1;
    }

    printf("Number of DOF : %ld\n", m.vtx_count);
    
    MinimalGraphSolver solver(m, test_f);

    solver.do_iterate_Newton(500,1e-05,1e-1,1e-6);
    solver.do_iterate_Picardi();
    

    return 0;

}