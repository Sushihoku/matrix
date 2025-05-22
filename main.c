#include "manip.h"
#include "matrix.h"
#include "oper.h"
#include <stdio.h>

int main() {

  // Testing the matrix exponent
  printf("Testing matrix exponential:\n");
  matrix *m = matrix_alloc(3, 3);
  *matrix_ptr(m, 0, 0) = 1.0;
  *matrix_ptr(m, 1, 1) = 2.0;
  *matrix_ptr(m, 2, 2) = -1.0;

  printf("Input matrix:\n");
  matrix_print(m);

  matrix *exp_m = matrix_exp(m, 1e-10);
  printf("Exponential of matrix:\n");
  matrix_print(exp_m);

  matrix_free(m);
  matrix_free(exp_m);

  printf("\nTesting Gauss method:\n");
  matrix *A = matrix_alloc(3, 3);
  *matrix_ptr(A, 0, 0) = 0;
  *matrix_ptr(A, 0, 1) = 0;
  *matrix_ptr(A, 0, 2) = 1.0;
  *matrix_ptr(A, 1, 0) = 0;
  *matrix_ptr(A, 1, 1) = 1.0;
  *matrix_ptr(A, 1, 2) = 0;
  *matrix_ptr(A, 2, 0) = 1.0;
  *matrix_ptr(A, 2, 1) = 0;
  *matrix_ptr(A, 2, 2) = 0;

  matrix *B = matrix_alloc(1, 3);
  *matrix_ptr(B, 0, 0) = 8.0;
  *matrix_ptr(B, 1, 0) = -11.0;
  *matrix_ptr(B, 2, 0) = -3.0;

  printf("Matrix A:\n");
  matrix_print(A);
  printf("\nRight part B:\n");
  matrix_print(B);

  matrix *X = matrix_solve_gauss(A, B);
  printf("\nSolution X:\n");
  matrix_print(X);

  // Checking the solution: calculate A*X - B
  if (X) {
    matrix *AX = matrix_alloc(1, 3);
    matrix_mul2(AX, A, X);
    matrix_sub(AX, B);
    printf("Residual A*X - B (should be near zero):\n");
    matrix_print(AX);
    printf("Residual norm: %f\n", matrix_norm(AX));
    matrix_free(AX);
  }

  matrix_free(A);
  matrix_free(B);
  matrix_free(X);

  return 0;
}
