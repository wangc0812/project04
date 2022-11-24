#ifndef  __MATRIX_H__
#define  __MATRIX_H__

#include <stdio.h>
#include "define.h"
#include <stdlib.h>

typedef struct Matrix {
    size_t row;
    size_t column;
    float *data;
} Matrix;

Matrix* createMatrix(size_t row, size_t column, size_t elenum, float* data);

void deleteMatrix(Matrix* mat);

void printMatrix(const Matrix* mat);

Matrix* copyMatrix(const Matrix* mat_src);

void array_sum(const float* a, const float* b, float* sum, size_t length);

void array_sub(const float* a, const float* b, float* sub, size_t length);

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

#endif