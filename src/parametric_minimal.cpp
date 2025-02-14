#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <functional>

#include "parametric_minimal.h"

#include "P1.h"
#include "tiny_blas.h"
#include "conjugate_gradient.h"
#include "minres.h"

#include "export_mesh.h"

#define HUGE 1e30

ParametricMinimalSolver::ParametricMinimalSolver(const Mesh &m, TVec3<std::function<double(const double)>> func, 
    TVec3<std::function<double(const double)>> dfunc, TVec3<std::function<double(const double)>> ddfunc)
    : m(m), N(m.vtx_count), N_b(m.boundary_count), u(N),A(N,N), J(N_b + 3,N_b + 3), b(N_b + 3), s(N_b), e_i(N_b, N), Se_i(N_b,N), inited(false),
    iterate(0), r(N,0.0), residual(iter_max, 0.0), p(N,0.0), Ap(N,0.0), converged(false){

    build_P1_CSRPattern(m, P);
    build_P1_stiffness_matrix(m, P, S);
    build_P1_stiffness_matrix(m, P, S0);

    for (size_t i = 0; i < N_b; i++){
        S0(m.boundary[i], m.boundary[i]) = HUGE;
        s[i] = 2.0*M_PI*i/N_b;
    }

    f = func;
    df = dfunc;
    ddf = ddfunc;

    for(size_t i = 0; i < N_b; i++){
        J(N_b,i) = 1.0;
        J(i,N_b) = 1.0;
    }

    J(N_b + 1,0) = -cos(s[1]) + 2*cos(s[0]) - cos(s[N_b - 1]);
    J(N_b + 2,0) = -sin(s[1]) + 2*sin(s[0]) - sin(s[N_b - 1]);
    J(0,N_b + 1) = -cos(s[1]) + 2*cos(s[0]) - cos(s[N_b - 1]);
    J(0,N_b + 2) = -sin(s[1]) + 2*sin(s[0]) - sin(s[N_b - 1]);

    for(size_t i = 1; i < N_b - 1; i++){
        J(N_b + 1,i) = cos(s[i + 1]) - 2*cos(s[i]) + cos(s[i - 1]);
        J(N_b + 2,i) = sin(s[i + 1]) - 2*sin(s[i]) + sin(s[i - 1]);
        J(i,N_b + 1) = cos(s[i + 1]) - 2*cos(s[i]) + cos(s[i - 1]);
        J(i,N_b + 2) = sin(s[i + 1]) - 2*sin(s[i]) + sin(s[i - 1]);
    }

    J(N_b + 1,N_b - 1) = -cos(s[0]) + 2*cos(s[N_b - 1]) - cos(s[N_b - 2]);
    J(N_b + 2,N_b - 1) = -sin(s[0]) + 2*sin(s[N_b - 1]) - sin(s[N_b - 2]);
    J(N_b - 1,N_b + 1) = -cos(s[0]) + 2*cos(s[N_b - 1]) - cos(s[N_b - 2]);
    J(N_b - 1,N_b + 2) = -sin(s[0]) + 2*sin(s[N_b - 1]) - sin(s[N_b - 2]);

    TArray<double> e_i_tmp(N,0.0), se_i(N, 0.0), r(N,0.0);
    e_i_tmp[m.boundary[0]] = -1.0;
    double errorCG = 0.0;
    for (size_t i = 0; i < N_b; i++) {
        S.mvp(e_i_tmp.data, se_i.data);
        conjugate_gradient_solve(S0, se_i.data, &e_i(i,0), r.data, p.data, Ap.data, &errorCG , toll*toll*toll, 10000, false);
        e_i_tmp[m.boundary[i]] = 0.0;
        if (i < N_b - 1)
            e_i_tmp[m.boundary[i + 1]] = -1.0;
        e_i(i,m.boundary[i]) += 1.0;
        S.mvp(&e_i(i,0), &Se_i(i,0));
    }

}

void ParametricMinimalSolver::clear_solution(){
    memset(u[0].data, 0.0, N * sizeof(double));
    memset(u[1].data, 0.0, N * sizeof(double));
    memset(u[2].data, 0.0, N * sizeof(double));
    converged = false;
    iterate = 0;
}

void ParametricMinimalSolver::do_iterate(size_t iter_max, double tol){
    
    clear_solution();
    TArray<double> gamma(N,0.0);
    TArray<double> Sgamma(N), S0gamma(N);
    TVec3<TArray<double>> Su(N);
    TArray<double> ds(N_b + 3);
    double alpha;
    double c = 1e-4;
    double rho = 0.5;
    TArray<double> s_tmp(N_b);
    bool flag;
    size_t iterCG;
    
    double error2 = 0.0, errorCG = 0.0, energy, energy_tmp = 0.0;

    for (size_t k = 0; k < 3; k++){
        for (size_t i = 0; i < N_b; i++)
            gamma[m.boundary[i]] = -f[k](s[i]);
        S.mvp(gamma.data, Sgamma.data);
        conjugate_gradient_solve(S0, Sgamma.data, u[k].data, r.data, p.data, Ap.data, &errorCG, toll, 10000, false);
        for(size_t i = 0; i < N_b; i++)
            u[k][m.boundary[i]] -= gamma[m.boundary[i]];
    }

    export_mesh_to_csv(m, u, "startParametric.csv");

    printf("Starting Parametric Minimal Solver.... \n");
    printf("%-8s %-15s %-8s %-15s %-15s\n", "Iter", "ErrorNewton", "IterCG", "ErrorCG", "Energy");

    while (iterate < iter_max){

        // reset gradient and Jacobian
       
        J.reset(N_b);
        memset(b.data, 0.0, N_b * sizeof(double));
        memset(ds.data, 0.0, (N_b + 3) * sizeof(double));

        // compute u

        for (size_t k = 0; k < 3; k++){
            for (size_t i = 0; i < N_b; i++)
                gamma[m.boundary[i]] = -f[k](s[i]);
            S.mvp(gamma.data, Sgamma.data);
            conjugate_gradient_solve(S0, Sgamma.data, u[k].data, r.data, p.data, Ap.data, &errorCG, toll, 10000, false);
            for(size_t i = 0; i < N_b; i++)
                u[k][m.boundary[i]] -= gamma[m.boundary[i]];
        }

        energy = 0.0;
        for(int k = 0; k < 3; k++){
            S.mvp(u[k].data, Su[k].data);
            energy += 0.5*blas_dot(u[k].data, Su[k].data, N);
        } 

        // compute rhs and Jacobian 
        for(int k = 0; k < 3; k++){
            for(size_t i = 0; i < N_b; i++){
                b[i] -= blas_dot(Su[k].data, &e_i(i,0), N)*df[k](s[i]);
                J(i,i) += blas_dot(Su[k].data, &e_i(i,0), N)*ddf[k](s[i]);
                for(size_t j = 0; j < N_b; j++)
                    J(i,j) += blas_dot(&Se_i(i,0), &e_i(j,0),N)*df[k](s[i])*df[k](s[j]);
            }
        }

        iterCG = minres_solve(J,b.data,ds.data,r.data,p.data,Ap.data, &errorCG, tol, 500, false);

        error2 = blas_dot(ds.data,ds.data, N_b);

        //Line search
        alpha = 1.0;
        flag = true;

        while(flag){
            for(size_t i = 0; i < N_b; i++)
                s_tmp[i] = s[i] + alpha*ds[i];
            energy_tmp = 0.0;
            for(int k = 0; k < 3; k++){
                for(size_t i = 0; i < N_b; i++)
                    gamma[m.boundary[i]] = -f[k](s_tmp[i]);
                S.mvp(gamma.data, Sgamma.data);
                conjugate_gradient_solve(S0, Sgamma.data, u[k].data, r.data, p.data, Ap.data, &errorCG, toll, 10000, false);
                for(size_t i = 0; i < N_b; i++)
                    u[k][m.boundary[i]] -= gamma[m.boundary[i]];
            }
            for(int k = 0; k < 3; k++){
                S.mvp(u[k].data, Su[k].data);
                energy_tmp += 0.5*blas_dot(u[k].data, Su[k].data, N);
            }
            if (energy_tmp < energy - c*alpha*error2 || alpha < 0.5)
                flag = false;
            else
                alpha *= rho;
        }

        
        for (size_t i = 0; i < N_b; i++)
            s[i] += alpha*ds[i];
        error2 = sqrt(error2);
        printf("%-8ld %-15g %-8ld %-15g %-15g\n", iterate, error2, iterCG, errorCG, energy);

        residual[iterate] = error2;

        if (error2 < tol)
            break;
        iterate++;
        
    }

    if (error2 <= tol) {
        converged = true;
        printf("Converged after %ld iterations.\n", iterate);
    }
    else {
        printf("Did not converge after %ld iterations.\n", iterate);
    }

    for (size_t k = 0; k < 3; k++){
        for (size_t i = 0; i < N_b; i++)
            gamma[m.boundary[i]] = -f[k](s[i]);
        S.mvp(gamma.data, Sgamma.data);
        conjugate_gradient_solve(S0, Sgamma.data, u[k].data, r.data, p.data, Ap.data, &errorCG, toll, 10000, false);
        for(size_t i = 0; i < N_b; i++)
            u[k][m.boundary[i]] -= gamma[m.boundary[i]];
    }
    export_mesh_to_csv(m, u, "solutionParametric.csv");
    
}