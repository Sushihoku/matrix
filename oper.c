#include "oper.h"
#include "matrix.h"

struct matrix {
  double *data;
  size_t w;
  size_t h;
};

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
