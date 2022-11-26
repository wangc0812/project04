#include "matrix.h"
#include "define.h"


Matrix* createMatrix(size_t row, size_t column, int elenum, float* data)
{
    //Generate Matrix Struct 
    //Please remember to free the memory due to dynamic allocate
    
    if (row <= 0 || column <= 0)
    {
        ERROR_INPUT_INPUTPARA;
        printf("ERROR: This error happened in 'createMatrix' \n");
        return NULL;
    }

    if (data == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'createMatrix' \n");
        return NULL;
    }
    int mat_size = row * column;
    if (mat_size == elenum)
    {
        Matrix* mat = MALLOC(1, Matrix);
        mat->row = row;
        mat->column = column;
        mat->data = MALLOC(mat_size, float);
        if (mat == NULL || mat->data == NULL)
        {
            free(mat);
            mat = NULL;

            free(mat->data);
            mat->data = NULL;

            ERROR_MEM_ALLOCATE;
            printf("ERROR: This error happened in 'createMatrix' \n");
            
            return NULL;
        }

        int i;
        for (i = 0; i < mat_size; i++) 
        {
            mat->data[i] = data[i];
        }

        return mat;
  
    }
    else
    {
        printf("ERROR: the number of data does not match the matrix size \n");
        printf("ERROR: This error happened in 'createMatrix' \n");
        
        return NULL;
    }
    
}

void deleteMatrix(Matrix* mat)
{
    if (mat == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'deleteMatrix()' \n");
        return;
    }

    FREE(mat -> data); 
    FREE(mat);
}

void printMatrix(const Matrix* mat) 
{
    //Print Matrix

    if (mat == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'printMatrix()' \n");
        return;
    }

    printf("---------------matrix print start---------------- \n");
    int i, j;
    for (i = 0; i < mat->row; i++) 
    {
        for (j = 0; j < mat->column; j++)
        {
            printf(PRECISION, mat->data[i * (mat->column) + j]);
        }
        printf("\n");
    }
    printf("---------------matrix print end---------------- \n");

    return;
}

Matrix* copyMatrix(const Matrix* mat_src) 
{
    //Copy Matrix
    if (mat_src == NULL)
    {
       ERROR_INPUT_POINTER;
       printf("ERROR: This error happened in 'copyMatrix()' \n");
       return NULL;
    }

    int elenum = mat_src->row * mat_src->column;
    Matrix* mat_copy = createMatrix(mat_src->row, mat_src->column, elenum, mat_src->data);
    
    return mat_copy;
}

void array_sum(const float* a, const float* b, float* sum, int length)
{
    if (a == NULL || b == NULL || sum == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'array_sum()' \n");
        return;
    }

    int  i;
    for (i = 0; i < length; i++)
    {
        sum[i] = a[i] + b[i];
    }

    return;
}

void array_sub(const float* a, const float* b, float* sub, int length)
{
    if (a == NULL || b == NULL || sub == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'array_sub()' \n");
        return;
    }

    int  i;
    for (i = 0; i < length; i++)
    {
        sub[i] = a[i] - b[i];
    }

    return;
}

Matrix* addMatrix(const Matrix* A, const Matrix* B)
{
    if (A == NULL || B == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'addMatrix()' \n");
        return NULL;
    }

    if (A->column != B->column || A->row != B->row)
    {
        ERROR_SIZE_MATCH;
        printf("ERROR: This error happened in 'addMatrix()' \n");
        return NULL;
    }

    int elenum = A->row * A->column;

    float sum[elenum];
    array_sum(A->data, B->data, sum, elenum);

    Matrix* mat_sum =  createMatrix(A->row, A->column, elenum, sum);

    return mat_sum;  
}

Matrix* subtractMatrix(const Matrix* A, const Matrix* B)
{
    if (A == NULL || B == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'subtractMatrix()' \n");
        return NULL;
    }

    if (A->column != B->column || A->row != B->row)
    {
        ERROR_SIZE_MATCH;
        printf("ERROR: This error happened in 'subtractMatrix()' \n");
        return NULL;
    }

    int elenum = A->row * A->column;

    float sub[elenum];
    array_sub(A->data, B->data, sub, elenum);

    Matrix* mat_sub = createMatrix(A->row, A->column, elenum, sub);

    return mat_sub;
}

Matrix* addScalar(const Matrix* A, const float b)
{
    if (A == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'addScalar()' \n");
        return NULL;
    }

    Matrix* C = copyMatrix(A);
    
    int  i;
    for (i = 0; i < (A->column * A->row); i++)
    {
        C->data[i] = C->data[i] + b;
    }

    return C;

}

Matrix* subScalar(const Matrix* A, const float b)
{
    if (A == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'subScalar()' \n");
        return NULL;
    }

    Matrix* C = addScalar(A, -1 * b);

    return C;

}

Matrix* multiplyScalar(const Matrix* A, const float b)
{
    if (A == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'multiplyScalar()' \n");
        return NULL;
    }

    Matrix* C = copyMatrix(A);
    
    int  i;
    for (i = 0; i < (A->column * A->row); i++)
    {
        C->data[i] = C->data[i] * b;
    }
    return C;
}

float maxelem(const Matrix* A)
{
    if (A == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'maxelem()' \n");
        return 1;
    }

    float max = A->data[0];
    int i;
    for ( i = 0; i < (A->column * A->row); i++)
	{
		if (max < A->data[i])
		{
			max = A->data[i];
		}
	}

    return max;
}

float minelem(const Matrix* A)
{
    if (A == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'minelem()' \n");
        return 1;
    }

    float min = A->data[0];
    int i;
    for ( i = 0; i < (A->column * A->row); i++)
	{
		if (min > A->data[i])
		{
			min = A->data[i];
		}
	}

    return min;
}

Matrix* mulMatrix(const Matrix* A, const Matrix* B)
{
    if (A == NULL || B == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'mulMatrix()' \n");
        return NULL;
    }

    if (A->column != B->row)
    {
        ERROR_SIZE_MATCH;
        printf("ERROR: This error happened in 'mulMatrix()' \n");
        return NULL;
    }

    int i, j, k, size;
    size = A->row * B->column;
    float c_data[size];
    int C_index = 0;
    float C_element = 0;

    for(i=0; i < A->row; i++)
    {
        for(j=0; j < B->column; j++)
        {
            for(k=0; k <A->column; k++)
            {
                C_element += A->data[(i * A->column) + k] * B->data[(k * B->column) + j];
            }
            c_data[C_index] = C_element;
            C_index += 1;
            C_element = 0;
        }
    }
    
    Matrix* C = createMatrix( A->row, B->column, size, c_data);

    return C;

}

// transposed matrix

Matrix* transpMatrix(const Matrix* A)
{
    if (A == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'transpMatrix()' \n");
        return NULL;
    }

    int i, j, size;
    size = A->row * A->column;
    float *B_data = MALLOC(size, float);
    memset(B_data, 0.0, size * sizeof(float));
    // float B_data[size];
    int B_index = 0;
    for(i=0; i<A->column; i++){
        for(j=0; j<A->row; j++){
            B_data[B_index] = A->data[j*A->column + i];
            B_index += 1;
        }
    }
    
    Matrix* B = createMatrix( A->column, A->row, size, B_data);

    return B;

}

// identity matrix
Matrix* identityMatrix(const size_t side)
{
    size_t i;
    size_t size = side * side;
    // float one = 1.0, zero=0.0;
    float *data = MALLOC(size, float);
    memset(data, 0.0, size * sizeof(float));
    // float data[size];

    // for (j = 0; j < size; j++)//c or cpp language
	// {
	// 	data[j] = zero;
	// }

    for(i=0; i<side; i++){
        data[i*side + i] = 1.0;
    };
    
    Matrix* mat = createMatrix(side, side, size, data);

    return mat;

}

Matrix* matmul_plain_row(const Matrix* A, const Matrix* B)
{
    if (A == NULL || B == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'matmul_SIMD()' \n");
        return NULL;
    }

    if (A->column != B->row)
    {
        ERROR_SIZE_MATCH;
        printf("ERROR: This error happened in 'matmul_SIMD()' \n");
        return NULL;
    }

    int i, j, k, size;
    size = A->row * B->column;
    float *C_data = MALLOC(size, float);
    memset(C_data, 0.0, size * sizeof(float));

    for(i=0; i < A->row; i++)
    {
        for(j=0; j < B->column; j++)
        {
            for(k=0; k <A->column; k++)
            {
               C_data[j + i * B->column] += A->data[(i * A->column) + k] * B->data[(k * B->column) + j];
            }
        }
    }
    
    Matrix* C = createMatrix( A->row, B->column, size, C_data);

    return C;

}

Matrix* matmul_plain_col(const Matrix* A, const Matrix* B)
{
    if (A == NULL || B == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'matmul_SIMD()' \n");
        return NULL;
    }

    if (A->column != B->row)
    {
        ERROR_SIZE_MATCH;
        printf("ERROR: This error happened in 'matmul_SIMD()' \n");
        return NULL;
    }

    int i, j, k, size;
    size = A->row * B->column;
    float *C_data = MALLOC(size, float);
    memset(C_data, 0.0, size * sizeof(float));

    for(i = 0; i < B->column; i++)
    {
        for(j = 0; j < A->row; j++)
        {
            for(k = 0; k <B->row; k++)
            {
                 C_data[i + j * B->column] += A->data[(j * A->column) + k] * B->data[i + k * B->column];
            }
        }
    }
    
    Matrix* C = createMatrix( A->row, B->column, size, C_data);

    return C;

}

Matrix* matmul_openmp(const Matrix* A, const Matrix* B)
{
    if (A == NULL || B == NULL)
    {
        ERROR_INPUT_POINTER;
        printf("ERROR: This error happened in 'matmul_SIMD()' \n");
        return NULL;
    }

    if (A->column != B->row)
    {
        ERROR_SIZE_MATCH;
        printf("ERROR: This error happened in 'matmul_SIMD()' \n");
        return NULL;
    }

    int i, j, k, size;
    size = A->row * B->column;
    float *C_data = MALLOC(size, float);
    memset(C_data, 0.0, size * sizeof(float));

    for(i=0; i < A->row; i++)
    {
        for(j=0; j < B->column; j++)
        {
            for(k=0; k <A->column; k++)
            {
               C_data[j + i * B->column] += A->data[(i * A->column) + k] * B->data[(k * B->column) + j];
            }
        }
    }
    
    Matrix* C = createMatrix( A->row, B->column, size, C_data);

    return C;
}