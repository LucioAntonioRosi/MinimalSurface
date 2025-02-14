#ifndef MINIMAL_GRAPH_SOLVER_H
#define MINIMAL_GRAPH_SOLVER_H

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <functional>

#include "minimal_graph.h"

#include "P1.h"
#include "tiny_blas.h"
#include "conjugate_gradient.h"
#include "minres.h"

#include "export_mesh.h"

#define HUGE 1e30

MinimalGraphSolver::MinimalGraphSolver(const Mesh &m, std::function<double(const Vec2d&)> func)
    : m(m), N(m.vtx_count), N_b(m.boundary_count), u(N), uold(N), du(N), q(m.tri_count), f(N),
        b(N, 0.0), r(N), p(N), Ap(N), inited(false), iterate_N(0), iterate_P(0), residual_N(iter_max,0.0), 
        residual_P(iter_max,0.0), converged(false){
    
    build_P1_CSRPattern(m, P); //
    
    for (size_t i = 0; i < N; i++) {
        f[i] = func(m.vertices[i]);
    }

}

void MinimalGraphSolver::clear_solution(bool Newton)
{
    memset(b.data, 0, N * sizeof(double));
    memset(du.data, 0, N * sizeof(double));
    for (size_t i = 0; i < N_b; i++) {
        b[m.boundary[i]] = f[m.boundary[i]];
    }
    memcpy(u.data, b.data, N * sizeof(double));
    memcpy(uold.data, b.data, N * sizeof(double));
    converged = false;
    if(Newton)
        iterate_N = 0;
    else
        iterate_P = 0;
}

double MinimalGraphSolver::compute_denominator(TArray<double> &den, const TArray<double> &u){
    uint32_t a, b, c;
    Vec2d A, B, C, AB, AC;
    double S_loc[6], area = 0.0, tri_area, ABAB, ACAC, ABAC, mult;
    for(int t = 0; t < m.tri_count; t++){
        a = m.triangles[t].x;
        b = m.triangles[t].y;
        c = m.triangles[t].z;
        A = m.vertices[a];
		B = m.vertices[b];
		C = m.vertices[c];
		AB = {(double)B[0] - (double)A[0],
			    (double)B[1] - (double)A[1]};
		AC = {(double)C[0] - (double)A[0],
			    (double)C[1] - (double)A[1]};

        ABAB = dot(AB,AB);
	    ACAC = dot(AC,AC);
	    ABAC = dot(AB, AC);
        tri_area = 0.5 * sqrt(ABAB * ACAC - ABAC * ABAC);
	    mult = 0.25 / (tri_area);
	    ABAB *= mult;
	    ACAC *= mult;
	    ABAC *= mult;

	    S_loc[0] = ACAC - 2 * ABAC + ABAB;
	    S_loc[1] = ACAC;
	    S_loc[2] = ABAB;
        S_loc[3] = ABAC - ACAC;
	    S_loc[4] = -ABAC;
	    S_loc[5] = ABAC - ABAB;

        den[t] = 1.0 / sqrt(1 + u[a]*u[a]*S_loc[0] + u[b]*u[b]*S_loc[1] + u[c]*u[c]*S_loc[2] +
                            2*(u[a]*u[b]*S_loc[3] + u[b]*u[c]*S_loc[4] + u[c]*u[a]*S_loc[5]));
        area += tri_area/den[t];
    }
    return area;
}



void MinimalGraphSolver::do_iterate_Newton(size_t max_iter, double tol, const double min_alpha, const double c, const double rho)
{

    clear_solution(true);
    TArray<double> u_tmp(N, 0.0);
    int iterCG;
    double error2 = 0.0, errorCG = 0.0, area;
    double alpha = 1.0, energy_tmp = 0.0;
    bool flag = true;

    memcpy(u.data,b.data ,N * sizeof(double));
    export_mesh_to_csv(m, u, "start.csv");

    area = compute_denominator(q, u);
    printf("Starting Newton solver.... \n");
    printf("%-10s %-15s %-15s %-15s %-15s\n", "Iter", "ErrorNewton", "IterCG", "ErrorCG", "Area");
    printf("%-10s %-15s %-15s %-15s %-15g \n", "-", "-", "-", "-", area);
    while (iterate_N < max_iter){
        
        build_P1_stiffness_matrix_NS(m, P, S_modified, q.data, u.data);
        build_P1_rhs_NS(m, q.data, u.data, b);

        for (size_t i = 0; i < N_b; i++){
            S_modified(m.boundary[i],m.boundary[i]) = HUGE;
            b[m.boundary[i]] = 0;           
        }

        iterCG = conjugate_gradient_solve(S_modified, b.data, du.data, r.data, p.data, Ap.data , &errorCG, tol, 10000, false);

        area = compute_denominator(q,u);
        error2 = blas_dot(du.data, du.data, N);

        flag = true;
        alpha = 1.0;

        while(flag){
            for (size_t i = 0; i < N; i++){
                u_tmp[i] = u[i] + alpha*du[i];
            }
            energy_tmp = compute_denominator(q,u_tmp);
            if (energy_tmp < area - c*alpha*error2 || alpha <= min_alpha){
                flag = false;
            }
            else{
                alpha *= rho;
            }
        }
        
        for (size_t i = 0; i < N; i++){
            u[i] += alpha*du[i];
        }

        
        error2 = sqrt(error2);
        
        printf("%-10ld %-15g %-15d %-15g %-15g\n", iterate_N, error2, iterCG, errorCG, area);

        residual_N[iterate_N] = error2;
        if (error2 < tol){
            break;
        }
        iterate_N++;
    }
    if (error2 <= tol) {
        converged = true;
        printf("Converged after %ld iterations.\n", iterate_N);
    }

    else {
        printf("Did not converge after %ld iterations.\n", iterate_N);
    }

    export_mesh_to_csv(m, u, "solutionNewton.csv");
}

void::MinimalGraphSolver::do_iterate_Picardi(size_t max_iter, double tol)
{   
    clear_solution(false);
    int iterCG;
    double error2 = 0.0, errorCG = 0.0, area;

    export_mesh_to_csv(m, u, "start.csv");
   
    for (size_t i = 0; i < N_b; i++){
        b[m.boundary[i]] *= HUGE;
    }

    area = compute_denominator(q,u);

    printf("Starting Picardi Solver.... \n");
    printf("%-10s %-15s %-15s %-15s %-15s\n", "Iter", "ErrorPicardi", "IterCG", "ErrorCG", "Area");
    printf("%-10s %-15s %-15s %-15s %-15g \n", "-", "-", "-", "-", area);
    while (iterate_P < max_iter){
        
        memcpy(uold.data, u.data, N * sizeof(double));

        build_P1_stiffness_matrix(m, P, S, true, q.data);

        for (size_t i = 0; i < N_b; i++){
            S(m.boundary[i],m.boundary[i]) = HUGE;           
        }

        iterCG = conjugate_gradient_solve(S, b.data, u.data, r.data, p.data, Ap.data , &errorCG, tol, 10000, false);

        error2 = 0.0;

        for(size_t i = 0; i < N; i++)
            error2 += (u[i] - uold[i])*(u[i] - uold[i]);
        
        error2 = sqrt(error2);

        area = compute_denominator(q,u);

        printf("%-10ld %-15g %-15d %-15g %-15g\n", iterate_P, error2, iterCG, errorCG, area);

        residual_P[iterate_P] = error2;

        if (error2 < tol){
            break;
        }
        iterate_P++;
    }
    
    if (error2 <= tol) {
        converged = true;
        printf("Converged after %ld iterations.\n", iterate_P);
    }

    else {
        printf("Did not converge after %ld iterations.\n", iterate_P);
    }
    export_mesh_to_csv(m, u, "solutionPicardi.csv");
}

#endif // MINIMAL_GRAPH_SOLVER_H