#include <math.h>
#include <algorithm>
#include <fstream>

#include "my_mesh.h"

void build_square_vertices(struct TVec2<double> *vert, size_t N, double a, double b){
    
    size_t V = N + 1;
    assert(V > 0);
    
    for (size_t i = 0; i < V; i++) {
        for (size_t j = 0; j < V; j++) {
            vert[V*i + j].x = (double)a*i/N;
            vert[V*i + j].y = (double)b*j/N;
        }
    }
    return;
}

void build_square_triangles(struct TVec3<size_t> *tri, size_t N){
    
    size_t V = N + 1, v, t = 0;
    
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            v = V*i + j;
            tri[t++] = {v, v + 1, v + 1 + V};
            tri[t++] = {v, v + 1 + V, v + V};
        }
    }
    return;
}

void square_boundary_edges(size_t *boundary, size_t N){
    
    size_t V = N + 1, k = 0;

    boundary[k++] = 0;
    boundary[k++] = N;
    boundary[k++] = V*N;
    boundary[k++] = V*N + N;
    
    for (size_t i = 1; i < N; i++) {
        boundary[k++] = i;
        boundary[k++] = V*i;
        boundary[k++] = V*i + N;
        boundary[k++] = V*N + i;
    }
    std::sort(boundary, boundary + k);

    return;
}

void build_square_mesh(struct Mesh *m, size_t N, double a, double b){
    
    size_t V = N + 1;
    m->boundary_count = 4*N;
    m->vtx_count  = V * V;
    m->vertices = (struct TVec2<double> *)malloc(m->vtx_count * sizeof(struct TVec2<double>));
    build_square_vertices(m->vertices, N, a, b);
    
    m->tri_count = 2 * N * N;
    m->triangles = (struct TVec3<size_t> *)malloc(m->tri_count * sizeof(struct TVec3<size_t>));
    build_square_triangles(m->triangles, N);
    
    m->boundary = (size_t *)malloc((4 * N) * sizeof(size_t));
    square_boundary_edges(m->boundary, N);
    
    return;
}


void build_disk_vertices(struct TVec2<double> *vert, size_t N, double R){
    
    assert(N > 0);
    size_t idx = 0;
    vert[idx++] = {0.0,0.0};
    for (size_t i = 1; i <= N ; i++) 
        for(size_t j = 0; j < 6*i; j++)
            vert[idx++] = {(double)R*i/N*cos(M_PI*j/(3*i)),(double)R*i/N*sin(M_PI*j/(3*i))};
    return;
}



void build_disk_triangles(struct TVec3<size_t> *tri, size_t N){
    
    size_t idx = 0, aux;
    tri[idx++] = {1, 2, 0};
    tri[idx++] = {2, 3, 0};
    tri[idx++] = {3, 4, 0};
    tri[idx++] = {4, 5, 0};
    tri[idx++] = {5, 6, 0};
    tri[idx++] = {6, 1, 0};

    for (size_t i = 2; i <= N; i++) {
        aux = 0;
        for (size_t j = 0; j < 6*i - 1; j++) {

            if (j % i != 0) {
                tri[idx++] = {2 + 3*(i-2)*(i-1) + aux,1 + 3*(i-1)*(i-2) + aux, 1 + 3*(i-1)*i + j};
                tri[idx++] = {2 + 3*(i-2)*(i-1) + aux,1 + 3*(i-1)*i + j, 2 + 3*(i-1)*i + j};
                ++aux;
            } else {
                tri[idx++] = {1 + 3*(i-1)*i + j, 2 + 3*(i-1)*i + j, 1 + 3*(i-2)*(i - 1) + aux};
            }
        }
        tri[idx++] = {1 + 3*(i-1)*i, 3*(i-1)*i + 6*i, 2 + aux + (i - 1)*(3*(i - 2) - 6)};
        tri[idx++] = {2 + aux + (i - 1)*(3*(i - 2) - 6), 1 + 3*(i-2)*(i - 1) + aux, 3*(i-1)*i + 6*i};
    }
    return;
}
 

void disk_boundary_edges(size_t *boundary, size_t N){
    
    for (size_t i = 0; i < 6*N; i++) 
        boundary[i] = 1 + 3*(N-1)*N + i;
    
    return;
}

void build_disk_mesh(struct Mesh *m, size_t N, double R){
    
    m->vtx_count = 1 + 3*N*(N+1);

    m->vertices = (struct TVec2<double> *) malloc(m->vtx_count * sizeof(struct TVec2<double>));
    build_disk_vertices(m->vertices, N, R);
    
    m->tri_count = 6*N*N;
    m->triangles = (struct TVec3<size_t> *)malloc(m->tri_count * sizeof(struct TVec3<size_t>));
    build_disk_triangles(m->triangles, N);

    m->boundary_count = 6*N;
    m->boundary = (size_t *)malloc(m->boundary_count * sizeof(size_t));
    disk_boundary_edges(m->boundary, N);
    
    return;
}
