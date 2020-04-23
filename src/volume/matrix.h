#include "global.h"

#ifndef MATRIX_H
#define MATRIX_H

typedef struct Matrix {
   FT* data; // size m * n
   // layout (row-wise):
   // a00, a01, ... a0n-1, [buffering]
   // ...
   // am-10, am-11, ... am-1n-1, [buffering] 
   int line; // size of one row + b, plus buffer for allignment
   int n; // dimensions
   int m; // constraints   
    
} Matrix;



Matrix* Matrix_new(int n, int m);


void Matrix_free(const void* o);

FT* Matrix_get_row(const Matrix* m, int i);
void Matrix_set(Matrix* m, int i, int x, FT a);
FT Matrix_get(const Matrix* m, int i, int x);
void Matrix_print(const void* o);
void Matrix_L_solve(const Matrix* o, FT* x, const FT* b);
void Matrix_invert_pdsym(const Matrix *In, Matrix *Out);

#endif
