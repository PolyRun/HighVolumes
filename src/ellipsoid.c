#include "beta_cut.h"



/**
 * \brief almost in-place cholesky factorization (from the book numerical recipes in c)
 * \param A a symmetric, positive definite matrix, the lower part of A will hold the subdiagonal entries of L
 * \param D will hold the diagonal of L
 * \param n the number of rows and cols of A
 **/
inline int cholesky(FT *A, FT *D, int n){
    for (int i = 1; i <= n; i++){
        for (int j = i; j <= n; j++){
            for (FT sum = A[i*n+j], int k = i-1; k >= 1; k--){
                sum -= A[i*n+k] * A[j*n+k];
            }
            if (i == j){
                // A is not positive definite (maybe due to rounding errors)
                if (sum <= 0){
                    return 1;
                }
                D[i] = sqrt(sum);
            }
            else {
                A[j*n + i] = sum/D[i];
            }
        }
    }

}




void Preprocess(Polytope *P, Polytope *Q, double *det){

    // MB: maybe implement this function as in PolyVest
    //checkHPs();

    int n = P->n;
    int m = P->m;

    
    double beta_r = 2*n;    
    double c3 = beta_r * beta_r;
    double c1 = (2 * n*n + (1-n/beta_r)*(1-n/beta_r)) * (1 - 1.0 / c3) / (2 * n*n - 2);
    double c2 = (1 - n / beta_r) / (n + 1);
    double c4 = 2 * c2 / (1 - 1.0 / beta_r);

    //init E(R2I, 0), T = R2I, ori = 0.
    
    FT R2;
    FT *ori;
    initEllipsoid(P, &R2, &ori);

    // initialize T to diag(R2)
    // T is out initial guess to the ellipsoid around Poly
    FT *T = (FT *) calloc(n*n, sizeof(FT));
    for (int i = 0; i < n; i++){
        T[i * n + i] = R2;
    }

    FT *distance = (FT *) calloc(m, sizeof(FT));
    FT *tm = (FT *) calloc(m, sizeof(FT));

    // update T until x^t*T*x <= 1
    //
    int counter = 0;
    while (++counter > 0){
        int i;
		
        //check if ori in polytope
        for(i = 0; i < m; i++) {
            FT sum = 0;
            for(int x=0; x < n; x++) {
                sum += ori[x] * Polytope_get_a(P, i, x);
            }
            // sum = (A * ori)_i > b_i (constraint not satisfied)
            if(sum > Polytope_get_b(P, i)) {
                // tm[i] = row_i(A)*T*row_i(A)^t
                tm[i] = 0;
                for (int j = 0; j < n; j++){
                    for (int k = 0; k < n; k++){
                        tm[i] += T[j*n + k] * Polytope_get_a(P, i, k);
                    }
                    tm[i] *= Polytope_get_a(P, i, j);
                }
                break;
            }
        }

        // check if small ellipsoid is contained in polytope
        if (i == m) {
            for (i = 0; i < m; i++){
                tm[i] = 0;
                for (int j = 0; j < n; j++){
                    FT tmi_tmp = 0;
                    for (int k = 0; k < n; k++){
                        tmi_tmp += T[j*n + k] * Polytope_get_a(P, i, k);
                    }
                    tm[i] += tmi_tmp * Polytope_get_a(P, i, j);
                }
                if (c3 * distance[i] * distance[i] - tm[i] < 0){
                    break;
                }
            }
        }
		
        //terminate if E satisfies the two criteria 
        if (i == m){
            break;
        }

        vec t = T * A.row(i).t() / sqrt(tm(i));
        FT *t = (FT *) malloc(n*sizeof(FT));
        for (int k = 0; k < n; k++){
            t[k] = 0;
            for (int j = 0; j < n; j++){
                t[k] += T[k*n+j] * Polytope_get_a(P, i, j);
            }
            t[k] /= sqrt(tm[i]);
        }
        for (int k = 0; k < n; k++){
            ori[k] -= t[k] * c2;
        }
        for (int k = 0; k < n; k++){
            for (int j = 0; j < n; j++){
                T[k*n + j] = c1 * (T[k*n + j] - c4 * t[k] * t[j]);  
            }
        }
    }


    /*
    if (!msg_off){ 
        cout << "R^2: " << R2 << endl << "Origin: " << endl;
    	ori.print();
    }
    */	

    
    //apply affine transformation in-place on Poly
    FT *D = (FT *) malloc(n*sizeof(FT));
    int err = cholesky(T, D, n);
    if (err > 0){
        printf("The input polytope is degenerated or non-existed and the volume is 0.\n");
        exit(1);		
    }

    // TODO from here down is not yet ported
    /*cout << Trans << endl;*/
    b = beta_r * (b - A * ori);
    A = A * Trans.t();

    
    printf("The number of iterations in shallow beta-cut: %d\n", counter);

    rowvec exp(n);
    exp.ones();
    for (int i = 0; i < n; i++){
        B[i] = b / A.col(i);
        Ai[i] = A / (A.col(i) * exp);
    }
	
    determinant = det(Trans) / pow(beta_r, n);
}
