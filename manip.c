#include "manip.h"
#include "oper.h"
#include <math.h>

struct matrix {
  double *data;
  size_t w;
  size_t h;
};

matrix *matrix_exp(const matrix *m, double eps) {
  if (!m || m->w != m->h)
    return NULL;

  size_t n = m->w;
  matrix *result = matrix_alloc_id(n, n);
  if (!result)
    return NULL;

  matrix *term = matrix_copy(m);
  if (!term) {
    matrix_free(result);
    return NULL;
  }

  int k = 1;
  while (matrix_norm(term) >= eps) {
    matrix *temp = matrix_copy(result);
    if (!temp) {
      matrix_free(result);
      matrix_free(term);
      return NULL;
    }

    if (matrix_add(temp, term) != 0) {
      matrix_free(temp);
      matrix_free(result);
      matrix_free(term);
      return NULL;
    }

    matrix_free(result);
    result = temp;

    // Calculate next term in the series: term = term * m / (k+1)
    matrix *next_term = matrix_alloc(n, n);
    if (!next_term) {
      matrix_free(result);
      matrix_free(term);
      return NULL;
    }

    if (matrix_mul2(next_term, term, m) != 0) {
      matrix_free(next_term);
      matrix_free(result);
      matrix_free(term);
      return NULL;
    }

    matrix_sdiv(next_term, ++k);

    matrix_free(term);
    term = next_term;
  }

  matrix_free(term);
  return result;
}

matrix *matrix_solve_gauss(const matrix *A, const matrix *B) {
  // Checking input parameters
  if (!A || !B || A->w != A->h || A->h != B->h || B->w != 1)
    return NULL;

  const size_t n = A->h;

  // Сreate copies of matrices A and B so as not to change the originals.
  matrix *ACopy = matrix_copy(A);
  matrix *BCopy = matrix_copy(B);
  if (!ACopy || !BCopy) {
    matrix_free(ACopy);
    matrix_free(BCopy);
    return NULL;
  }

  // Direct run of Gauss's method
  for (size_t k = 0; k < n; ++k) {
    // Finding the row with the maximum element in column k
    size_t max_row = k;
    double max_val = fabs(*matrix_ptr(ACopy, k, k));

    for (size_t i = k + 1; i < n; ++i) {
      double val = fabs(*matrix_ptr(ACopy, i, k));
      if (val > max_val) {
        max_val = val;
        max_row = i;
      }
    }

    // Rearranging string
    if (max_row != k) {
      matrix_swap_rows(ACopy, k, max_row);

      double tmp = *matrix_ptr(BCopy, k, 0);
      *matrix_ptr(BCopy, k, 0) = *matrix_ptr(BCopy, max_row, 0);
      *matrix_ptr(BCopy, max_row, 0) = tmp;
    }

    // Test for degeneracy
    if (fabs(*matrix_ptr(ACopy, k, k)) < 1e-12) {
      matrix_free(ACopy);
      matrix_free(BCopy);
      return NULL;
    }

    //  Excluding column k in bottom rows
    for (size_t i = k + 1; i < n; ++i) {
      double factor = *matrix_ptr(ACopy, i, k) / *matrix_ptr(ACopy, k, k);
      *matrix_ptr(BCopy, i, 0) -= factor * *matrix_ptr(BCopy, k, 0);

      for (size_t j = k; j < n; ++j) {
        *matrix_ptr(ACopy, i, j) -= factor * *matrix_ptr(ACopy, k, j);
      }
    }
  }

  // Reverse Gaussian Method
  matrix *X = matrix_alloc(1, n);
  if (!X) {
    matrix_free(ACopy);
    matrix_free(BCopy);
    return NULL;
  }

  for (size_t i = n; i-- > 0;) {
    double sum = 0.0;
    for (size_t j = i + 1; j < n; ++j) {
      sum += *matrix_ptr(ACopy, i, j) * *matrix_ptr(X, j, 0);
    }
    *matrix_ptr(X, i, 0) =
        (*matrix_ptr(BCopy, i, 0) - sum) / *matrix_ptr(ACopy, i, i);
  }

  matrix_free(ACopy);
  matrix_free(BCopy);
  return X;
}
