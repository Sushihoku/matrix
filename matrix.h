#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

// Incomplete matrix type (encapsulation)
struct matrix;
typedef struct matrix matrix;

// Matrix creation and destruction
matrix *matrix_alloc(size_t w, size_t h);
matrix *matrix_copy(const matrix *m);
void matrix_free(matrix *m);

// Fast element access
double *matrix_ptr(matrix *m, size_t i, size_t j);
const double *matrix_cptr(const matrix *m, size_t i, size_t j);

// Matrix initialization
void matrix_set_zero(matrix *m);
void matrix_set_id(matrix *m);
matrix *matrix_alloc_zero(size_t w, size_t h);
matrix *matrix_alloc_id(size_t w, size_t h);

// Assignment
int matrix_assign(matrix *m1, const matrix *m2);

// Matrix manipulations
void matrix_transpose(matrix *m);
void matrix_swap_rows(matrix *m, size_t i1, size_t i2);
void matrix_swap_cols(matrix *m, size_t j1, size_t j2);
void matrix_mul_row(matrix *m, size_t i, double d);
void matrix_add_rows(matrix *m, size_t i1, size_t i2);
double matrix_norm(const matrix *m);

// Input/output
void matrix_print(const matrix *m);
matrix *matrix_input(size_t w, size_t h);

#endif // MATRIX_Hendif
