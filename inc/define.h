#ifndef  __DEFINE_H__
#define  __DEFINE_H__

// typedef float DATA_TYPE;
// typedef int INT;
// typedef void VOID;

#define PRECISION "%.5lf\t"

#define MALLOC(n, type) \
         ((type *)malloc( (n) * sizeof(type)))

#define FREE(p)             \
         if(p!=NULL)         \
         {                   \
               free(p);       \
               p = NULL;     \
         }
#define ERROR_INPUT_POINTER printf("ERROR: input data pointer error \n")
#define ERROR_INPUT_INPUTPARA printf("ERROR:input paramter error \n")
#define ERROR_MEM_ALLOCATE printf("ERROR:failed to allocate memeory \n")
#define ERROR_SIZE_MATCH printf("ERROR:matrix size does not match \n")
#endif