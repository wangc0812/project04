#ifndef  __MATRIX_H__
#define  __MATRIX_H__

#include <stdio.h>
#include "define.h"
#include <stdlib.h>
#include <xmmintrin.h>
#include <string.h>
#include <omp.h>

typedef struct Matrix {
    size_t row;
    size_t column;
    float *data;
} Matrix;

Matrix* createMatrix(size_t row, size_t column, int elenum, float* data);

void deleteMatrix(Matrix* mat);

void printMatrix(const Matrix* mat);

Matrix* copyMatrix(const Matrix* mat_src);

void array_sum(const float* a, const float* b, float* sum, int length);

void array_sub(const float* a, const float* b, float* sub, int length);

Matrix* addMatrix(const Matrix* A, const Matrix* B);

Matrix* subtractMatrix(const Matrix* A, const Matrix* B);

Matrix* addScalar(const Matrix* A, const float b);

Matrix* subScalar(const Matrix* A, const float b);

Matrix* multiplyScalar(const Matrix* A, const float b);

float maxelem(const Matrix* A);

float minelem(const Matrix* A);

Matrix* mulMatrix(const Matrix* A, const Matrix* B);

Matrix* transpMatrix(const Matrix* A);

Matrix* identityMatrix(const size_t side);

Matrix* matmul_plain_row(const Matrix* A, const Matrix* B);

Matrix* matmul_plain_col(const Matrix* A, const Matrix* B);

Matrix* matmul_SIMD(const Matrix* A, const Matrix* B);

#endif