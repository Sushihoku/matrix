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
