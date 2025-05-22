#ifndef MANIP_H
#define MANIP_H
#include "matrix.h"

// Special functions
matrix *matrix_exp(const matrix *m, double eps);
matrix *matrix_solve_gauss(const matrix *A, const matrix *B);

#endif
