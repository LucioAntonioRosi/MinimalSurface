#include <stdio.h>
#include <string.h>

#ifdef USE_OPENMP
	#include <omp.h>
#endif

#include "P1.h"
#include "adjacency.h"
#include "my_mesh.h"
#include "sparse_matrix.h"
#include "vec2.h"

static bool find(uint32_t x, uint32_t *start, size_t count)
{
	for (size_t i = 0; i < count; ++i) {
		if (start[i] == x)
			return true;
	}
	return false;
}

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P)
{
	size_t vtx_count = m.vtx_count;
	size_t tri_count = m.tri_count;

	P.row_start.resize(vtx_count + 1);

	VTAdjacency adj(m);

	/* Upper bound on the number of edges a->b with a <= b */
	size_t max_nnz = 3 * tri_count + vtx_count;
	P.col.resize(max_nnz);

	/* Fill P.row_start and P.col (not yet ordered) */
	size_t nnz = 0;
	for (size_t a = 0; a < vtx_count; ++a) {
		P.row_start[a] = nnz;
		uint32_t *start = &P.col[nnz];
		size_t nnz_loc = 0;
		uint32_t kstart = adj.offset[a];
		uint32_t kstop = kstart + adj.degree[a];
		for (size_t k = kstart; k < kstop; ++k) {
			uint32_t b = adj.vtri[k].next;
			uint32_t c = adj.vtri[k].prev;
			if (b < a && !find(b, start, nnz_loc)) {
				P.col[nnz++] = b;
				nnz_loc++;
			}
			if (c < a && !find(c, start, nnz_loc)) {
				P.col[nnz++] = c;
				nnz_loc++;
			}
		}
		P.col[nnz++] = a;
	}
	P.row_start[vtx_count] = nnz;
	P.col.resize(nnz);
	P.col.shrink_to_fit();

	/* Reorder each "line" of P.col in increasing column order */
	for (size_t a = 0; a < vtx_count; ++a) {
		uint32_t *to_sort = &P.col[P.row_start[a]];
		size_t count = P.row_start[a + 1] - P.row_start[a];
		/* Insertion sort */
		for (size_t k = 1; k < count; ++k) {
			size_t j = k - 1;
			while (j && (to_sort[j] > to_sort[j + 1])) {
				uint32_t tmp = to_sort[j];
				to_sort[j] = to_sort[j + 1];
				to_sort[j + 1] = tmp;
				j--;
			}
		}
	}
}

static void stiffness(const Vec2d &AB, const Vec2d &AC, double *__restrict S);
static void stiffness_NS(const Vec2d &AB, const Vec2d &AC, double *__restrict S, const double *u, const double den);

void inline stiffness(const Vec2d &AB, const Vec2d &AC, double *__restrict S)
{
	double ABAB = dot(AB,AB);
	double ACAC = dot(AC,AC);
	double ABAC = dot(AB, AC);
	double mult = 0.5 / sqrt(ABAB * ACAC - ABAC * ABAC);
	ABAB *= mult;
	ACAC *= mult;
	ABAC *= mult;

	S[0] = ACAC - 2 * ABAC + ABAB;
	S[1] = ACAC;
	S[2] = ABAB;
	S[3] = ABAC - ACAC;
	/* Note the chosen order : (B,C)-> 4 and (C,A) -> 5 */
	S[4] = -ABAC;
	S[5] = ABAC - ABAB;

}

static void stiffness_NS(const Vec2d &AB, const Vec2d &AC, double *__restrict S, const double *u_loc, const double den){
	double S_start[6];
	stiffness(AB, AC, S_start);
	double GraduGradPhi[3];

	GraduGradPhi[0] = u_loc[0]*S_start[0] + u_loc[1]*S_start[3] + u_loc[2]*S_start[5];
	GraduGradPhi[1] = u_loc[0]*S_start[3] + u_loc[1]*S_start[1] + u_loc[2]*S_start[4];
	GraduGradPhi[2] = u_loc[0]*S_start[5] + u_loc[1]*S_start[4] + u_loc[2]*S_start[2];

	S[0] = S_start[0] - den*den*(GraduGradPhi[0]*GraduGradPhi[0]);
	S[1] = S_start[1] - den*den*(GraduGradPhi[1]*GraduGradPhi[1]);
	S[2] = S_start[2] - den*den*(GraduGradPhi[2]*GraduGradPhi[2]);
	S[3] = S_start[3] - den*den*(GraduGradPhi[0]*GraduGradPhi[1]);
	S[4] = S_start[4] - den*den*(GraduGradPhi[1]*GraduGradPhi[2]);
	S[5] = S_start[5] - den*den*(GraduGradPhi[2]*GraduGradPhi[0]);
}

void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S, bool modified, const double *den) 
{
	size_t vtx_count = m.vtx_count;
	size_t tri_count = m.tri_count;
	assert(P.row_start.size == vtx_count + 1);

	S.symmetric = true;
	S.rows = S.cols = vtx_count;
	S.nnz = P.col.size;
	S.row_start = P.row_start.data;
	S.col = P.col.data;
	S.data.resize(S.nnz);
	for (size_t i = 0; i < S.nnz; ++i) {
		S.data[i] = 0.0;
	}
	size_t a,b,c;
	Vec2d A, B, C, AB, AC;
	double Sloc[6];
	for (size_t t = 0; t < tri_count; ++t) {
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
		
		stiffness(AB, AC, Sloc);
		if (modified) {
			Sloc[0] *= den[t];
			Sloc[1] *= den[t];
			Sloc[2] *= den[t];
			Sloc[3] *= den[t];
			Sloc[4] *= den[t];
			Sloc[5] *= den[t];
		}

		S(a, a) += Sloc[0];
		S(b, b) += Sloc[1];
		S(c, c) += Sloc[2];
		S(a > b ? a : b, a > b ? b : a) += Sloc[3];
		S(b > c ? b : c, b > c ? c : b) += Sloc[4];
		S(c > a ? c : a, c > a ? a : c) += Sloc[5];
	}
}

void build_P1_stiffness_matrix_NS(const Mesh &m, const CSRPattern &P, CSRMatrix &S, const double *den, const double *u) 
{
	size_t vtx_count = m.vtx_count;
	size_t tri_count = m.tri_count;
	assert(P.row_start.size == vtx_count + 1);

	S.symmetric = true;
	S.rows = S.cols = vtx_count;
	S.nnz = P.col.size;
	S.row_start = P.row_start.data;
	S.col = P.col.data;
	S.data.resize(S.nnz);
	for (size_t i = 0; i < S.nnz; ++i) {
		S.data[i] = 0.0;
	}

	size_t a,b,c;
	Vec2d A, B, C, AB, AC;
	double Sloc[6], u_loc[3];

	for (size_t t = 0; t < tri_count; ++t) {
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
	
		u_loc[0] = u[a];
		u_loc[1] = u[b];
		u_loc[2] = u[c];
		stiffness_NS(AB, AC, Sloc, u_loc, den[t]);
		
		S(a, a) += Sloc[0]*den[t];
		S(b, b) += Sloc[1]*den[t];
		S(c, c) += Sloc[2]*den[t];
		S(a > b ? a : b, a > b ? b : a) += Sloc[3]*den[t];
		S(b > c ? b : c, b > c ? c : b) += Sloc[4]*den[t];
		S(c > a ? c : a, c > a ? a : c) += Sloc[5]*den[t];
	}
}

void build_P1_rhs_NS(const Mesh &m, const double *den, const double *u, TArray<double> &rhs){
	
	memset(rhs.data, 0, m.vtx_count * sizeof(double));
	double S_loc[6];
	size_t a,b,c;
	Vec2d A, B, C, AB, AC;
	for(size_t t = 0; t < m.tri_count; t++){
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

		stiffness(AB, AC, S_loc);
		rhs[a] -= den[t]*(u[a]*S_loc[0] +u[b]*S_loc[3] + u[c]*S_loc[5]);
		rhs[b] -= den[t]*(u[a]*S_loc[3] +u[b]*S_loc[1] + u[c]*S_loc[4]);
		rhs[c] -= den[t]*(u[a]*S_loc[5] +u[b]*S_loc[4] + u[c]*S_loc[2]);
	}

}