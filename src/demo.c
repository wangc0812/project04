#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "define.h"

int demo()
{
    printf("test \n");
    const int elenum = 6;
    float data[6] = {5,2,8,4,2,5};
    const int row = 2;
    const int col = 3;

    // generate matrix
    printf("**************************** \n");
    printf("create a Matrix \n");
    Matrix* mat1 = createMatrix(row, col, elenum, data);

    // print matrix
    printMatrix(mat1);

    // copy matrix
    printf("**************************** \n");
    printf("copy a Matrix \n");
    Matrix* mat2 = copyMatrix(mat1);
    printMatrix(mat2);  // to confirm the copy function

    printf("**************************** \n");
    printf("add two Matrices \n");
    Matrix* mat3 = addMatrix(mat1, mat2); // matrix add
    printMatrix(mat3);  // to confirm the add function
    printf("**************************** \n");
    printf("subtrucy two Matrices \n");
    Matrix* mat4 = subtractMatrix(mat3, mat2); // matrix sub
    printMatrix(mat4);  // to confirm the sub function

    // matrix add scalar
    printf("**************************** \n");
    printf("matrix add scalar(6.6) \n");
    float b = 6.6;
    Matrix* mat5 = addScalar(mat1, b);
    printMatrix(mat5);  // to confirm the add scalar function

    printf("**************************** \n");
    printf("matrix subtract scalar(6.6) \n");
    Matrix* mat6 = subScalar(mat1, b);
    printMatrix(mat6);  // to confirm the add scalar function

    printf("**************************** \n");
    printf("matrix mutiple scalar(2.0) \n");
    float c = 2.0;
    Matrix* mat7 = multiplyScalar(mat1, c);
    printMatrix(mat7);  // to confirm the add scalar function

    printf("**************************** \n");
    printf("find the minimal and maximal values of a matrix \n");
    float  max = maxelem(mat1);
    float  min = minelem(mat1);
    printf("max value = %f, min value = %f\n", max , min );

    printf("**************************** \n");
    printf("multiply two matrices with unmatched size. \n");
    Matrix* mat8 = mulMatrix(mat1, mat2);
    printMatrix(mat8);

    printf("**************************** \n");
    printf("transpose matrices. \n");
    Matrix* mat9 = transpMatrix(mat1);
    printMatrix(mat9);

    printf("**************************** \n");
    printf("multiply two matrices after transpose a Matrix. \n");
    Matrix* mat10 = mulMatrix(mat1, mat9);
    printMatrix(mat10);

    printf("**************************** \n");
    printf("create an identity Matrix. \n");
    int side = 8;
    Matrix* mat11 = identityMatrix(side);
    printMatrix(mat11);


    // free memory
    deleteMatrix(mat1);
    deleteMatrix(mat2);
    deleteMatrix(mat3);
    deleteMatrix(mat4);
    deleteMatrix(mat5);
    deleteMatrix(mat6);
    deleteMatrix(mat7);
    deleteMatrix(mat8);
    deleteMatrix(mat9);
    deleteMatrix(mat10);
    deleteMatrix(mat11);
    
    return 0;
}