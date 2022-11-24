#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "define.h"
int main()
{
    printf("test \n");
    const int elenum = 6;
    float data[6] = {5,2,8,4,2,5};
    const int row = 2;
    const int col = 3;

    // generate matrix
    printf("create a Matrix \n");
    Matrix* mat1 = createMatrix(row, col, elenum, data);

    // print matrix
    printMatrix(mat1);

    Matrix* mat2 = transpMatrix(mat1); 
    printMatrix(mat2);  // to confirm the add function


    Matrix* mat3 = matmul_plain(mat1, mat2);
    printMatrix(mat3);


    
    // free memory
    deleteMatrix(mat1);
    deleteMatrix(mat2);
    deleteMatrix(mat3);

    
    return 0;
}