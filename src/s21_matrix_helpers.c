#include "s21_matrix_helpers.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "s21_matrix.h"

int allocate_memory(int rows, int columns, matrix_t* matrix) {
  int error = OK;
  double** tmp = (double**)malloc(rows * columns * sizeof(double) +
                                  rows * sizeof(double*));

  if (tmp == NULL) {
    error = INCORRECT_MATRIX;
  } else {
    double* pointer = (double*)(tmp + rows);
    for (int i = 0; i < rows; i++) {
      *(tmp + i) = pointer + columns * i;
    }
  }
  matrix->matrix = tmp;
  return error;
}

void free_memory(matrix_t* matrix) {
  if (matrix->matrix != NULL) {
    free(matrix->matrix);
    matrix->matrix = NULL;
  }
}

void init_valid_matrix(matrix_t* matrix) {
  for (int i = 0; i < matrix->rows; i++) {
    for (int j = 0; j < matrix->columns; j++) {
      matrix->matrix[i][j] = 0.;
    }
  }
}

int check_matrix(matrix_t* matrix) {
  int error = SUCCESS;

  if (NULL == matrix || matrix->matrix == NULL || matrix->columns < 1 ||
      matrix->rows < 1) {
    error = FAILURE;
  }
  return error;
}

void swap_strings(matrix_t* matrix, int str1, int str2) {
  for (int i = 0; i < matrix->columns; i++) {
    double tmp = 0;
    tmp = matrix->matrix[str1][i];
    matrix->matrix[str1][i] = matrix->matrix[str2][i];
    matrix->matrix[str2][i] = tmp;
  }
}

void fill_minor_matrix(matrix_t* matrix, int row, int column,
                       matrix_t* minor_matrix) {
  for (int i = 0, new_row = 0; i < matrix->rows; i++) {
    if (row != i) {
      for (int j = 0, new_column = 0; j < matrix->rows; j++) {
        if (column != j) {
          minor_matrix->matrix[new_row][new_column++] = matrix->matrix[i][j];
        }
      }
      new_row++;
    }
  }
}

void write_random_matrix(matrix_t* matrix) {
  srand(time(NULL));

  for (int i = 0; i < matrix->rows; i++) {
    for (int j = 0; j < matrix->columns; j++) {
      matrix->matrix[i][j] =
          (double)(rand()) / RAND_MAX * pow(-1., rand() % 2) * (rand() % 100);
    }
  }
}

void write_iteratively_matrix(matrix_t* matrix, double first, double step) {
  for (int i = 0; i < matrix->rows; i++) {
    for (int j = 0; j < matrix->columns; j++) {
      matrix->matrix[i][j] = first;
      first += step;
    }
  }
}