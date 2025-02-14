#include <iostream>
#include <string.h>
#include <functional>

#include "parametric_minimal.h"
#include "my_mesh.h"
#include "tiny_blas.h"
#include "P1.h"

int main(int argc, char **argv){

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <number of vertices subcircles for the disk>" << std::endl;
        std::cerr << "Mesh type: disk only" << std::endl;
        return 1;
    }

    int N = atoi(argv[1]);
    double R = 1.5;

    if (N <= 0) {
        std::cerr << "Number of vertices must be a positive integer." << std::endl;
        return 1;
    }

    Mesh m;
    
    build_disk_mesh(&m, N);

    printf("Number of DOF : %ld\n", m.vtx_count);

    TVec3<std::function<double(const double)>> func, dfunc, ddfunc;
    //func[0] = [R](const double x) { return (1 + 0.25*cos(3*x))*cos(2*x); };
    //func[1] = [R](const double x) { return (1 + 0.25*cos(3*x))*sin(2*x); };
    //func[2] = [R](const double x) { return 0.25*sin(3*x); };
    func[0] = [R](const double x) { return cos(x); };
    func[1] = [R](const double x) { return sin(x); };
    //func[2] = [R](const double x) { return sin(6*M_PI*cos(x)) + sin(2*M_PI*sin(x)); };
    func[2] = [R](const double x) {return 5*exp(cos(x)) + 1/(1 + sin(x)*sin(x)) + 2*cos(x)*sin(x) + pow(5*cos(x),2);};
    dfunc[0] = [R, func](const double x) { return (func[0](x + 0.001) - func[0](x))/0.001; };
    dfunc[1] = [R, func](const double x) { return (func[1](x + 0.001) - func[1](x))/0.001; };
    dfunc[2] = [R, func](const double x) { return (func[2](x + 0.001) - func[2](x))/0.001; };
    ddfunc[0] = [R, dfunc](const double x) { return (dfunc[0](x + 0.001) - dfunc[0](x))/0.001; };
    ddfunc[1] = [R, dfunc](const double x) { return (dfunc[1](x + 0.001) - dfunc[1](x))/0.001; };
    ddfunc[2] = [R, dfunc](const double x) { return (dfunc[2](x + 0.001) - dfunc[2](x))/0.001; };
    //func[0] = [R](const double x) { return R*R*R*cos(3*x) + 4*R*R*R*R*R*cos(5*x); };
    //func[1] = [R](const double x) { return R*R*R*sin(3*x) - 4*R*R*R*R*R*sin(5*x); };
    //func[2] = [R](const double x) { return -sqrt(15)*R*R*R*R*sin(4*x); };
    //dfunc[0] = [R](const double x) { return -3*R*R*R*sin(3*x) - 20*R*R*R*R*R*sin(5*x); };
    //dfunc[1] = [R](const double x) { return 3*R*R*R*cos(3*x) - 20*R*R*R*R*R*cos(5*x); };
    //dfunc[2] = [R](const double x) { return -4*sqrt(15)*R*R*R*R*cos(4*x); };
    //ddfunc[0] = [R](const double x) { return -9*R*R*R*cos(3*x) - 100*R*R*R*R*R*cos(5*x); };
    //ddfunc[1] = [R](const double x) { return -9*R*R*R*sin(3*x) + 100*R*R*R*R*R*sin(5*x); };
    //ddfunc[2] = [R](const double x) { return 16*sqrt(15)*R*R*R*R*sin(4*x); };

    //func[0] = [R](const double x) { return R*cos(x) + R*R*R/3*cos(3*x); };
    //func[1] = [R](const double x) { return R*sin(x) - R*R*R/3*sin(3*x); };
    //func[2] = [R](const double x) { return R*R*cos(2*x); };
    //dfunc[0] = [R](const double x) { return -R*sin(x) - R*R*R*sin(3*x); };
    //dfunc[1] = [R](const double x) { return R*cos(x) - R*R*R*cos(3*x); };
    //dfunc[2] = [R](const double x) { return -2*R*R*sin(2*x); };
    //ddfunc[0] = [R](const double x) { return -R*cos(x) - 3*R*R*R*cos(3*x); };
    //ddfunc[1] = [R](const double x) { return -R*sin(x) + 3*R*R*R*sin(3*x); };
    //ddfunc[2] = [R](const double x) { return -4*R*R*cos(2*x); };

    ParametricMinimalSolver solver(m, func, dfunc, ddfunc);

    solver.do_iterate();
    
    return 0;

}