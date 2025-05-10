#include "matrix.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct matrix {
  double *data;
  size_t w;
  size_t h;
};

matrix *matrix_alloc(size_t w, size_t h) {
  matrix *m = malloc(sizeof(matrix));
  if (!m)
    return NULL;

  m->data = calloc(
      w * h,
      sizeof(
          double)); // Используем calloc вместо malloc для инициализации нулями
  if (!m->data) {
    free(m);
    return NULL;
  }

  m->w = w;
  m->h = h;
  return m;
}

void matrix_free(matrix *m) {
  if (m) {
    free(m->data);
    free(m);
  }
}

matrix *matrix_copy(const matrix *m) {
  if (!m)
    return NULL;

  matrix *copy = matrix_alloc(m->w, m->h);
  if (!copy)
    return NULL;

  memcpy(copy->data, m->data, m->w * m->h * sizeof(double));
  return copy;
}

double *matrix_ptr(matrix *m, size_t i, size_t j) {
  return &m->data[i * m->w + j];
}

const double *matrix_cptr(const matrix *m, size_t i, size_t j) {
  return &m->data[i * m->w + j];
}

void matrix_set_zero(matrix *m) {
  if (!m)
    return;
  memset(m->data, 0, m->w * m->h * sizeof(double));
}

void matrix_set_id(matrix *m) {
  if (!m)
    return;
  matrix_set_zero(m);
  size_t min_dim = m->w < m->h ? m->w : m->h;
  for (size_t i = 0; i < min_dim; ++i) {
    *matrix_ptr(m, i, i) = 1.0;
  }
}

matrix *matrix_alloc_zero(size_t w, size_t h) {
  matrix *m = matrix_alloc(w, h);
  if (m)
    matrix_set_zero(m);
  return m;
}

matrix *matrix_alloc_id(size_t w, size_t h) {
  matrix *m = matrix_alloc(w, h);
  if (m)
    matrix_set_id(m);
  return m;
}

int matrix_assign(matrix *m1, const matrix *m2) {
  if (!m1 || !m2 || m1->w != m2->w || m1->h != m2->h)
    return -1;

  memcpy(m1->data, m2->data, m1->w * m1->h * sizeof(double));
  return 0;
}

void matrix_transpose(matrix *m) {
  if (!m)
    return;

  // In-place transposition only for square matrices
  if (m->w == m->h) {
    for (size_t i = 0; i < m->h; ++i) {
      for (size_t j = i + 1; j < m->w; ++j) {
        double tmp = *matrix_ptr(m, i, j);
        *matrix_ptr(m, i, j) = *matrix_ptr(m, j, i);
        *matrix_ptr(m, j, i) = tmp;
      }
    }
  } else {
    // For non-square matrices, create a temporary copy
    matrix *temp = matrix_copy(m);
    if (!temp)
      return;

    free(m->data);
    m->data = malloc(m->w * m->h * sizeof(double));
    if (!m->data) {
      m->data = temp->data;
      temp->data = NULL;
      matrix_free(temp);
      return;
    }

    size_t old_w = m->w;
    size_t old_h = m->h;
    m->w = old_h;
    m->h = old_w;

    for (size_t i = 0; i < m->h; ++i) {
      for (size_t j = 0; j < m->w; ++j) {
        *matrix_ptr(m, i, j) = *matrix_cptr(temp, j, i);
      }
    }

    matrix_free(temp);
  }
}

void matrix_swap_rows(matrix *m, size_t i1, size_t i2) {
  if (!m || i1 >= m->h || i2 >= m->h)
    return;

  for (size_t j = 0; j < m->w; ++j) {
    double tmp = *matrix_ptr(m, i1, j);
    *matrix_ptr(m, i1, j) = *matrix_ptr(m, i2, j);
    *matrix_ptr(m, i2, j) = tmp;
  }
}

void matrix_swap_cols(matrix *m, size_t j1, size_t j2) {
  if (!m || j1 >= m->w || j2 >= m->w)
    return;

  for (size_t i = 0; i < m->h; ++i) {
    double tmp = *matrix_ptr(m, i, j1);
    *matrix_ptr(m, i, j1) = *matrix_ptr(m, i, j2);
    *matrix_ptr(m, i, j2) = tmp;
  }
}

void matrix_mul_row(matrix *m, size_t i, double d) {
  if (!m || i >= m->h)
    return;

  for (size_t j = 0; j < m->w; ++j) {
    *matrix_ptr(m, i, j) *= d;
  }
}

void matrix_add_rows(matrix *m, size_t i1, size_t i2) {
  if (!m || i1 >= m->h || i2 >= m->h)
    return;

  for (size_t j = 0; j < m->w; ++j) {
    *matrix_ptr(m, i1, j) += *matrix_ptr(m, i2, j);
  }
}

double matrix_norm(const matrix *m) {
  if (!m || m->w == 0 || m->h == 0)
    return 0.0;

  double max_sum = 0.0;
  for (size_t i = 0; i < m->h; ++i) {
    double row_sum = 0.0;
    for (size_t j = 0; j < m->w; ++j) {
      row_sum += fabs(*matrix_cptr(m, i, j));
    }
    if (row_sum > max_sum) {
      max_sum = row_sum;
    }
  }
  return max_sum;
}

void matrix_print(const matrix *m) {
  if (!m) {
    printf("NULL matrix\n");
    return;
  }

  for (size_t i = 0; i < m->h; ++i) {
    for (size_t j = 0; j < m->w; ++j) {
      printf("%8.4f ", *matrix_cptr(m, i, j));
    }
    printf("\n");
  }
}

matrix *matrix_input(size_t w, size_t h) {
  matrix *m = matrix_alloc(w, h);
  if (!m)
    return NULL;

  printf("Enter matrix %zux%zu:\n", h, w);
  for (size_t i = 0; i < h; ++i) {
    for (size_t j = 0; j < w; ++j) {
      scanf("%lf", matrix_ptr(m, i, j));
    }
  }

  return m;
}
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
int matrix_add(matrix *m1, const matrix *m2) {
  if (!m1 || !m2 || m1->w != m2->w || m1->h != m2->h)
    return -1;

  for (size_t i = 0; i < m1->h; ++i) {
    for (size_t j = 0; j < m1->w; ++j) {
      *matrix_ptr(m1, i, j) += *matrix_cptr(m2, i, j);
    }
  }
  return 0;
}

int matrix_sub(matrix *m1, const matrix *m2) {
  if (!m1 || !m2 || m1->w != m2->w || m1->h != m2->h)
    return -1;

  for (size_t i = 0; i < m1->h; ++i) {
    for (size_t j = 0; j < m1->w; ++j) {
      *matrix_ptr(m1, i, j) -= *matrix_cptr(m2, i, j);
    }
  }
  return 0;
}

void matrix_smul(matrix *m, double d) {
  if (!m)
    return;

  for (size_t i = 0; i < m->h; ++i) {
    for (size_t j = 0; j < m->w; ++j) {
      *matrix_ptr(m, i, j) *= d;
    }
  }
}

void matrix_sdiv(matrix *m, double d) {
  if (!m || d == 0.0)
    return;
  matrix_smul(m, 1.0 / d);
}

int matrix_add2(matrix *m, const matrix *m1, const matrix *m2) {
  if (!m || !m1 || !m2 || m1->w != m2->w || m1->h != m2->h || m->w != m1->w ||
      m->h != m1->h)
    return -1;

  for (size_t i = 0; i < m->h; ++i) {
    for (size_t j = 0; j < m->w; ++j) {
      *matrix_ptr(m, i, j) = *matrix_cptr(m1, i, j) + *matrix_cptr(m2, i, j);
    }
  }
  return 0;
}

int matrix_sub2(matrix *m, const matrix *m1, const matrix *m2) {
  if (!m || !m1 || !m2 || m1->w != m2->w || m1->h != m2->h || m->w != m1->w ||
      m->h != m1->h)
    return -1;

  for (size_t i = 0; i < m->h; ++i) {
    for (size_t j = 0; j < m->w; ++j) {
      *matrix_ptr(m, i, j) = *matrix_cptr(m1, i, j) - *matrix_cptr(m2, i, j);
    }
  }
  return 0;
}

int matrix_smul2(matrix *m, const matrix *m1, double d) {
  if (!m || !m1 || m->w != m1->w || m->h != m1->h)
    return -1;

  for (size_t i = 0; i < m->h; ++i) {
    for (size_t j = 0; j < m->w; ++j) {
      *matrix_ptr(m, i, j) = *matrix_cptr(m1, i, j) * d;
    }
  }
  return 0;
}

int matrix_sdiv2(matrix *m, const matrix *m1, double d) {
  if (d == 0.0)
    return -1;
  return matrix_smul2(m, m1, 1.0 / d);
}

int matrix_mul(matrix *m1, const matrix *m2) {
  if (!m1 || !m2 || m1->w != m2->h)
    return -1;

  matrix *temp = matrix_alloc(m2->w, m1->h);
  if (!temp)
    return -1;

  for (size_t i = 0; i < m1->h; ++i) {
    for (size_t j = 0; j < m2->w; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < m1->w; ++k) {
        sum += *matrix_cptr(m1, i, k) * *matrix_cptr(m2, k, j);
      }
      *matrix_ptr(temp, i, j) = sum;
    }
  }

  int result = matrix_assign(m1, temp);
  matrix_free(temp);
  return result;
}

int matrix_mul2(matrix *m, const matrix *m1, const matrix *m2) {
  if (!m || !m1 || !m2 || m1->w != m2->h || m->w != m2->w || m->h != m1->h)
    return -1;

  // If m is the same as m1 or m2, we need to use a temporary matrix
  if (m == m1 || m == m2) {
    matrix *temp = matrix_alloc(m2->w, m1->h);
    if (!temp)
      return -1;

    int result = matrix_mul2(temp, m1, m2);
    if (result == 0) {
      result = matrix_assign(m, temp);
    }
    matrix_free(temp);
    return result;
  }

  for (size_t i = 0; i < m1->h; ++i) {
    for (size_t j = 0; j < m2->w; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < m1->w; ++k) {
        sum += *matrix_cptr(m1, i, k) * *matrix_cptr(m2, k, j);
      }
      *matrix_ptr(m, i, j) = sum;
    }
  }
  return 0;
}
