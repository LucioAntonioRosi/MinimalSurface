#pragma once

#include <functional>

#include "sparse_matrix.h"
#include "matrix.h"
#include "my_mesh.h"

#include "vec2.h"
#include "vec3.h"
#include "array.h"

struct ParametricMinimalSolver {
    ParametricMinimalSolver(const Mesh &m, TVec3<std::function<double(const double)>> func,
     TVec3<std::function<double(const double)>> dfunc, TVec3<std::function<double(const double)>> ddfunc);
    const Mesh &m;
    size_t N;   
    size_t N_b; 

    TVec3<TArray<double>> u;
    CSRPattern P;
    CSRMatrix S, S0; 
    Matrix A, J;

    TVec3<std::function<double(const double)>> f, df, ddf; 
    TArray<double> b, s;
    Matrix e_i, Se_i;

    bool inited; 
    size_t iterate; 
    TArray<double> r, residual, p, Ap; 
    bool converged; 

    static const size_t iter_max = 500;
    static constexpr double toll = 1e-5;

    void clear_solution();
    
    void do_iterate(size_t max_iter = iter_max, double tol = toll);
    void do_iterate2(size_t max_iter = iter_max, double tol = toll);
};