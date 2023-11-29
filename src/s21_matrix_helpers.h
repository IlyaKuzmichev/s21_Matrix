#pragma once

#include "s21_matrix.h"

int allocate_memory(int rows, int columns, matrix_t* matrix);
void free_memory(matrix_t* matrix);
void init_valid_matrix(matrix_t* matrix);
int check_matrix(matrix_t* matrix);
void swap_strings(matrix_t* matrix, int str1, int str2);
void fill_minor_matrix(matrix_t* matrix, int row, int column,
                       matrix_t* minor_matrix);
void write_random_matrix(matrix_t* matrix);
void write_iteratively_matrix(matrix_t* matrix, double first, double step);
