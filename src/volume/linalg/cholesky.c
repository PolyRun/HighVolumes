#include "cholesky.h"


/**
 * \brief cholesky factorization
 * \param A a symmetric, positive definite matrix
 * \param Trans will be lower diagonal matrix holding L
 * \param n the number of rows and cols of A
 **/
int cholesky(FT *A, FT *Trans, int n);


inline int cholesky(FT *A, FT *Trans, int n){
    
    // Decomposing a matrix into Lower Triangular 
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j <= i; j++) { 
            FT sum = 0; 
  
            if (j == i) { 
                for (int k = 0; k < j; k++) {
                    sum += Trans[j*n+k] * Trans[j*n+k];
                }
                if (A[j*n+j] < sum){
                    // matrix not pd
                    return 1;
                }
                Trans[j*n+j] = sqrt(A[j*n+j] - sum); 
            }
            else { 
                for (int k = 0; k < j; k++) {
                    sum += (Trans[i*n+k] * Trans[j*n+k]);
                }
                Trans[i*n+j] = (A[i*n+j] - sum) / Trans[j*n+j]; 
            } 
        } 
    }
    return 0;
}


int cholesky_matrix(const Matrix *M, Matrix *L){
    int n = M->n;
    
    // Decomposing the matrix M into Lower Triangular 
    for (int i = 0; i < n; i++) {
        FT *Mi = Matrix_get_row(M, i);
        for (int j = 0; j <= i; j++) { 
            FT sum = 0;

            // L_{j,j} = sqrt(T_{j,j} - sum_{k = 0}^{j-1} L_{j, k}^2)
            if (j == i) {
                FT *Lj = Matrix_get_row(L, j);
                sum = dotProduct(Lj, Lj, j);
                if (Mi[j] < sum){
                    // matrix not pd
                    return 1;
                }
                Matrix_set(L, j, j, sqrt(Mi[j] - sum));
            }
            // for j < i
            // L_{i,j} = 1/L_{j,j} * (T_{i,j} - sum_{k = 0}^{j-1} l_{i,k} * l_{j, k})
            else {
                FT *Li = Matrix_get_row(L, i);
                FT *Lj = Matrix_get_row(L, j);
                FT sum = dotProduct(Li, Lj, j);
                Matrix_set(L, i, j, 1/Lj[j] * (Mi[j] - sum));
            } 
        } 
    }
    return 0;
}


int cholesky_ellipsoid(const Ellipsoid *E, Matrix *L){
    int n = E->n;
    
    // Decomposing the matrix T of E into Lower Triangular 
    for (int i = 0; i < n; i++) {
        FT *Ti = Ellipsoid_get_Ti(E, i);
        for (int j = 0; j <= i; j++) { 
            FT sum = 0;

            // L_{j,j} = sqrt(T_{j,j} - sum_{k = 0}^{j-1} L_{j, k}^2)
            if (j == i) {
                FT *Lj = Matrix_get_row(L, j);
                sum = dotProduct(Lj, Lj, j);
                if (Ti[j] <= sum){
                    Matrix_print(L);
                    printf("%f <= %f\n", Ti[j], sum);
                    return 1;
                }
                Matrix_set(L, j, j, sqrt(Ti[j] - sum));
            }
            // for j < i
            // L_{i,j} = 1/L_{j,j} * (T_{i,j} - sum_{k = 0}^{j-1} l_{i,k} * l_{j, k})
            else {
                FT *Li = Matrix_get_row(L, i);
                FT *Lj = Matrix_get_row(L, j);
                FT sum = dotProduct(Li, Lj, j);
                Matrix_set(L, i, j, 1.0/Lj[j] * (Ti[j] - sum));
            } 
        } 
    }
    return 0;
}



