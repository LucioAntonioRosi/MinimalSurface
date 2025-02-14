#pragma once

#include "sparse_matrix.h"
#include "my_mesh.h"
#include "array.h"

void build_P1_CSRPattern(const Mesh &m, CSRPattern &P);
void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P,
			       CSRMatrix &S, bool modified = false, const double *den = nullptr);
void build_P1_stiffness_matrix_NS(const Mesh &m, const CSRPattern &P,
				   CSRMatrix &S, const double *den, const double *u);
void build_P1_rhs_NS(const Mesh &m, const double *den, const double *u,
				   TArray<double> &rhs);