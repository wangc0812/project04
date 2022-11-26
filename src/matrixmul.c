// #include <stdio.h>
// #include <stdlib.h>
// #include "matrix.h"
// #include "define.h"
// #include <xmmintrin.h>


// int main()
// {
    
//     // printf("test \n");
//     // const int elenum = 6;
//     // float data[6] = {5,2,8,4,2,5};
//     // const int row = 2;
//     // const int col = 3;

//     // // generate matrix
//     // printf("create a Matrix \n");
//     // Matrix* mat1 = createMatrix(row, col, elenum, data);

//     // // print matrix
//     // printMatrix(mat1);

//     // Matrix* mat2 = transpMatrix(mat1); 
//     // printMatrix(mat2);  // to confirm the add function


//     // Matrix* mat3 = matmul_plain_row(mat1, mat2);
//     // printMatrix(mat3);

//     // Matrix* mat4 = matmul_plain_col(mat1, mat2);
//     // printMatrix(mat4);
    
//     // // free memory
//     // deleteMatrix(mat1);
//     // deleteMatrix(mat2);
//     // deleteMatrix(mat3);
//     // deleteMatrix(mat4);

//     // test simd
//     // float a[4] = { 1,2,3,4 };
//     // float b[4] = { 5,6,7,8 };
//     // float res[4];
//     // // __mm_mm，表示这是一个64位/128位的指令，_mm256和_mm512则表示是256位或是512位的指令
//     // //_loadu，表示unaligen的load指令，不带u后缀的为aligen版本
//     // //u，表示 unaligned，内存未对齐。如果是a，表示 aligned，内存已对齐
//     // //p，表示 packed，打包数据，会对128位所有数据执行操作。如果是s，则表示 scalar，标量数据，仅对128位内第一个数执行操作
//     // //s，表示 single precision floating point，将数据视为32位单精度浮点数，一组4个
//     // //如果是d，表示 double precision floating point，将数据视为64位双精度浮点，一组两个。
//     // //_ps，同上面汇编指令，还可以是_pd，_ss，_sd

//     // __m128 A = _mm_loadu_ps(a);
//     // __m128 B = _mm_loadu_ps(b);
//     // __m128 RES = _mm_mul_ps(A, B);
//     // //__m128变量本身是占128位内存的结构体，不开编译器优化的时候会频繁发生回写内存再读取，性能很差
//     // _mm_storeu_ps(res, RES);


    
//     return 0;
// }

#include <stdio.h>
#include <omp.h>
int main(int argc, char **argv)
{
    int nthreads, thread_id;
    printf("I am the main thread.\n");
#pragma omp parallel private(nthreads, thread_id)
    {
        nthreads = omp_get_num_threads();
        thread_id = omp_get_thread_num();
        printf("Hello. I am thread %d out of a team of %d\n", thread_id, nthreads);
    }
    printf("Here I am, back to the main thread.\n");
    return 0;
}

