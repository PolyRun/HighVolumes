#include "matrix.h"

Matrix* Matrix_new(int n, int m) {
   Matrix* o = (Matrix*) malloc(sizeof(Matrix));
   o->n = n;
   o->m = m;
   o->line = ceil_cache(n,sizeof(FT)); // make sure next is also 32 alligned
   o->data = (FT*)(aligned_alloc(32, o->line*m*sizeof(FT))); // align this to 32
   for(int i=0;i<(o->line*m);i++) {o->data[i]=0;}
   return o;
}


void Matrix_free(const void* o) {
   Matrix* p = (Matrix*)o;
   free(p->data);
   free(p);
}

inline FT* Matrix_get_row(const Matrix* m, int i) {
   return &(m->data[i * (m->line)]);
}

inline void Matrix_set(Matrix* m, int i, int x, FT a) {
   m->data[i * (m->line) + x] = a;
}
inline FT Matrix_get(const Matrix* m, int i, int x) {
   return m->data[i * (m->line) + x];
}

void Matrix_print(const void* o) {
   const Matrix* p = (Matrix*)o;
   printf("Matrix: n=%d, m=%d\n",p->n,p->m);
   for(int i=0; i<p->m; i++) {
      for(int j=0; j<p->n; j++) {
         printf(" %.3f", Matrix_get(p,i,j));
      }
      printf("\n");
   }
}

void Matrix_L_solve(const Matrix* o, FT* x, const FT* b) {
   const Matrix* L = (Matrix*)o;
   const int n = L->n;

   for(int i=0;i<n;i++) {
      FT sum = 0;
      FT* Li = Matrix_get_row(L, i);
      for(int j=0;j<i;j++) {
         sum += x[j] * Li[j];
      }
      //assert(abs(Li[i]) > 1e-8 && "numerically unstable");
      x[i] = (b[i] - sum) / Li[i];
   }
}



void Matrix_invert_pdsym(const Matrix *In, Matrix *Out){
    const int n = In->n;

    Matrix* L = Matrix_new(n,n);
    Matrix* Linvt = Matrix_new(n,n);
    FT* b = (FT*)(aligned_alloc(32, n*sizeof(FT))); // align this to 32
    int err = cholesky_matrix(In, L);
    assert(err==0 && "no cholesky errors");
   
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {b[j]=(i==j);}// unit vec
        FT* x = Matrix_get_row(Linvt, i);
        Matrix_L_solve(L, x, b);
    }
    //Matrix_print(Linvt);
   
    // A = T.inverse()
    // A = (LLt).inverse()
    // A = Lt.inverse() * L.inverse()
    for(int i=0;i<n;i++) {
        FT* Ai = Matrix_get_row(Out,i);
        for(int j=0;j<n;j++) {
            FT* a = Matrix_get_row(Linvt, i);
            FT* b = Matrix_get_row(Linvt, j);
            FT dot = dotProduct(a,b, n);
            Ai[j] = dot;
        }
    }

    free(b);
    Matrix_free(Linvt);
    Matrix_free(L);
}

void Matrix_rotate(Matrix *m, const int i, const int j, const FT angle) {
    if(i==j) {return;}
    const int n = m->n;
    FT s = sin(angle);
    FT c = cos(angle);
    for(int k=0;k<n;k++) {
        FT vi = Matrix_get(m, i, k);
        FT vj = Matrix_get(m, j, k);
	Matrix_set(m, i, k, c*vi - s*vj);
	Matrix_set(m, j, k, s*vi + c*vj);
    }
}

void Matrix_MVM(const Matrix *M, const FT* x, FT* y) {
   const int n=M->n;
   const int m=M->m;
   for(int i=0;i<m;i++) {
      FT yi = 0;
      for(int j=0;j<n;j++) {
         yi += Matrix_get(M,i,j) * x[j];
      }
      y[i] = yi;
      //printf("MVM %d %f\n",i,yi);
   }
}
