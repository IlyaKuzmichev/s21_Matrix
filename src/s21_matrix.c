#include "s21_matrix.h"

#include <math.h>
#include <stdlib.h>

#include "s21_matrix_helpers.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = OK;

  if (NULL == result || rows < 1 || columns < 1) {
    error = INCORRECT_MATRIX;
  } else {
    error = allocate_memory(rows, columns, result);
    result->rows = rows;
    result->columns = columns;
  }
  if (OK == error) {
    init_valid_matrix(result);
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (NULL != A) {
    free_memory(A);
    A->columns = A->rows = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int error = SUCCESS;

  if (!(check_matrix(A) && check_matrix(B) && A->rows == B->rows &&
        A->columns == B->columns)) {
    error = FAILURE;
  }
  for (int i = 0; i < A->rows && SUCCESS == error; i++) {
    for (int j = 0; j < A->columns && SUCCESS == error; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= S21_EPS) {
        error = FAILURE;
      }
    }
  }
  return error;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;

  if (!(check_matrix(A) && check_matrix(B) && result)) {
    error = INCORRECT_MATRIX;
  }
  if (OK == error && (A->rows != B->rows || A->columns != B->columns)) {
    error = CALCULATION_ERROR;
  }
  if (OK == error) {
    error = s21_create_matrix(A->rows, A->columns, result);
  }
  for (int i = 0; i < A->rows && OK == error; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;

  if (!(check_matrix(A) && check_matrix(B) && result)) {
    error = INCORRECT_MATRIX;
  }
  if (OK == error && (A->rows != B->rows || A->columns != B->columns)) {
    error = CALCULATION_ERROR;
  }
  if (OK == error) {
    error = s21_create_matrix(A->rows, A->columns, result);
  }
  for (int i = 0; i < A->rows && OK == error; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }
  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = OK;
  if (!(check_matrix(A) && result)) {
    error = INCORRECT_MATRIX;
  }
  if (OK == error) {
    error = s21_create_matrix(A->rows, A->columns, result);
  }
  for (int i = 0; i < A->rows && OK == error; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
    }
  }
  return error;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;

  if (!(check_matrix(A) && check_matrix(B) && result)) {
    error = INCORRECT_MATRIX;
  }
  if (OK == error && A->columns != B->rows) {
    error = CALCULATION_ERROR;
  }
  if (OK == error) {
    error = s21_create_matrix(A->rows, B->columns, result);
  }
  for (int i = 0; i < A->rows && OK == error; i++) {
    for (int j = 0; j < B->columns; j++) {
      for (int k = 0; k < A->columns; k++) {
        result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
      }
    }
  }
  return error;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = OK;

  if (!(check_matrix(A) && result)) {
    error = INCORRECT_MATRIX;
  }
  if (OK == error) {
    error = s21_create_matrix(A->columns, A->rows, result);
  }
  for (int i = 0; i < A->rows && OK == error; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[j][i] = A->matrix[i][j];
    }
  }
  return error;
}

int s21_determinant(matrix_t *A, double *result) {
  int error = OK;

  if (!(check_matrix(A) && result)) {
    error = INCORRECT_MATRIX;
  }
  if (OK == error && A->rows != A->columns) {
    error = CALCULATION_ERROR;
  }
  if (OK == error && 1 == A->rows) {
    *result = A->matrix[0][0];
  } else if (OK == error) {
    matrix_t temp_matrix = {0};
    *result = 1.;
    error = s21_create_matrix(A->rows, A->columns, &temp_matrix);
    for (int i = 0; i < A->rows && OK == error; i++) {
      for (int j = 0; j < A->columns; j++) {
        temp_matrix.matrix[i][j] = A->matrix[i][j];
      }
    }
    for (int i = 0;
         i < temp_matrix.rows && OK == error && fabs(*result) > S21_EPS; i++) {
      int swap_string = i;
      while (swap_string < A->rows - 1 &&
             fabs(temp_matrix.matrix[swap_string][i]) < S21_EPS) {
        swap_string++;
      }
      if (swap_string != i) {
        swap_strings(&temp_matrix, swap_string, i);
        *result *= -1;
      }
      if (fabs(temp_matrix.matrix[i][i]) >= S21_EPS) {
        *result *= temp_matrix.matrix[i][i];
      } else {
        *result = 0.;
      }
      for (int j = i + 1; j < temp_matrix.rows && fabs(*result) > S21_EPS;
           j++) {
        double multiplier = temp_matrix.matrix[j][i] / temp_matrix.matrix[i][i];
        for (int k = i; k < temp_matrix.rows; k++) {
          temp_matrix.matrix[j][k] -= multiplier * temp_matrix.matrix[i][k];
        }
      }
    }
    s21_remove_matrix(&temp_matrix);
  }
  return error;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = OK;

  if (!(check_matrix(A) && result)) {
    error = INCORRECT_MATRIX;
  }
  if (OK == error &&
      (A->rows != A->columns || 1 == A->rows || 1 == A->columns)) {
    error = CALCULATION_ERROR;
  }
  if (OK == error) {
    error = s21_create_matrix(A->columns, A->rows, result);
  }
  if (OK == error && A->rows > 1) {
    matrix_t minor_matrix = {0};
    error = s21_create_matrix(A->rows - 1, A->rows - 1, &minor_matrix);
    for (int i = 0; i < A->rows && OK == error; i++) {
      for (int j = 0; j < A->rows; j++) {
        fill_minor_matrix(A, i, j, &minor_matrix);
        double determinant = 0;
        s21_determinant(&minor_matrix, &determinant);
        result->matrix[i][j] = determinant * pow(-1., (double)(i + j));
      }
    }
    s21_remove_matrix(&minor_matrix);
  }
  return error;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = OK;

  if (!(check_matrix(A) && result)) {
    error = INCORRECT_MATRIX;
  } else {
    *result = (matrix_t){0};
  }
  if (OK == error && A->rows != A->columns) {
    error = CALCULATION_ERROR;
  }
  if (OK == error) {
    double determinant = 0.;
    error = s21_determinant(A, &determinant);
    if (fabs(determinant) < S21_EPS) {
      error = CALCULATION_ERROR;
    }
    if (OK == error) {
      if (1 == A->rows) {
        error = s21_create_matrix(1, 1, result);
        if (OK == error) {
          result->matrix[0][0] = 1 / determinant;
        }
      } else {
        matrix_t complement_matrix = {0};
        error = s21_calc_complements(A, &complement_matrix);
        if (OK == error) {
          matrix_t transponent_matrix = {0};
          error = s21_transpose(&complement_matrix, &transponent_matrix);
          if (OK == error) {
            error = s21_mult_number(&transponent_matrix,
                                    (double)(1. / determinant), result);
          }
          s21_remove_matrix(&transponent_matrix);
        }
        s21_remove_matrix(&complement_matrix);
      }
    }
  }
  return error;
}