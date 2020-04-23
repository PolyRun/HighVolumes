//#include "beta_cut.h"
#include "preprocess.h"
#include <assert.h>

#define FLOATWIDTH 15


// ------------ Find real implementation in cholesky.h, included through volume.h
///**
// * \brief cholesky factorization
// * \param A a symmetric, positive definite matrix
// * \param Trans will be lower diagonal matrix holding L
// * \param n the number of rows and cols of A
// **/
//int cholesky(FT *A, FT *Trans, int n);
//
//
//inline int cholesky(FT *A, FT *Trans, int n){
//    
//    // Decomposing a matrix into Lower Triangular
//    
//    //------------------------------
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j <= i; j++) { 
//            FT sum = 0; 
//
//            //------------------------------
//            if (j == i) { 
//                for (int k = 0; k < j; k++) {
//                    sum += Trans[j*n+k] * Trans[j*n+k];
//                }
//                Trans[j*n+j] = sqrt(A[j*n+j] - sum); 
//            }
//            // (i+1, i, 0, 1)
//            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//            
//            //------------------------------
//            else { 
//                for (int k = 0; k < j; k++) {
//                    sum += (Trans[i*n+k] * Trans[j*n+k]);
//                }
//                Trans[i*n+j] = (A[i*n+j] - sum) / Trans[j*n+j]; 
//            }
//            // (j+1, j, 1, 0)
//            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//        }
//    }
//    // sum_{i=0}^{n} (i+1, i, 0, 1) +
//    // sum_{i=1}^{n-1} sum_{j=0}^{i-1} (j+1, j, 1, 0)
//    // = ((n^3 + 6n^2 + 5n)/6,
//    //    (n^3 + 3n^2 - 4n)/6,
//    //    (n^2 + n)/2,
//    //    n)
//    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//    return 0;
//}



void preprocess(Polytope *P, Polytope **Q, FT *det){

    // MB: maybe implement this function as in PolyVest
    //checkHPs();

    int n = P->n;
    int m = P->m;

    *Q = Polytope_new(n, m);
    
    double beta_r = 2*n; 
    double c1 = (2 * pow(n, 2) + pow(1 - n / beta_r, 2)) * (1 - 1.0 / pow(beta_r, 2)) / (2 * pow(n, 2) - 2);
    double c2 = (1 - n / beta_r) / (n + 1);
    double c3 = beta_r * beta_r;
    double c4 = 2 * c2 / (1 - 1.0 / beta_r);

#ifdef DEBUG
    printf("c1: %0.*f\n", FLOATWIDTH, c1);
    printf("c2: %0.*f\n", FLOATWIDTH, c2);
    printf("c3: %0.*f\n", FLOATWIDTH, c3);
    printf("c4: %0.*f\n", FLOATWIDTH, c4);
#endif
    
    //init E(R2I, 0), T = R2I, ori = 0.
    
    FT R2;
    FT *ori;
    init_ellipsoid(P, &R2, &ori);


    // initialize T to diag(R2)
    // T is out initial guess to the ellipsoid around Poly
    FT *T = (FT *) calloc(n*n, sizeof(FT));
    for (int i = 0; i < n; i++){
        T[i * n + i] = R2;
    }

#ifdef DEBUG_MSG
    printf("--------------- HIGHVOLUMES\n");
    printf("First ellipsoid approximation\n");
    printf("\nT:\n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            printf("%0.*f ",FLOATWIDTH, T[i*n + j]);
        }
        printf("\n");
    }
    printf("\ncenter:\n");
    for (int i = 0; i < n; i++){
        printf("%0.*f ", FLOATWIDTH, ori[i]);
    }
    printf("\n");
#endif

    FT *distance = (FT *) calloc(m, sizeof(FT));
    FT *tm = (FT *) calloc(m, sizeof(FT));
    FT *t = (FT *) malloc(n*sizeof(FT));

    // update T until x^t*T*x <= 1

    int counter = 0;
    while (++counter > 0){
        int i;

        

        // distance = b - A * ori
#ifdef DEBUG
        //printf("ROUND %d\n", counter);
#endif
        for (i = 0; i < m; i++){
            FT sum = 0;
            for(int x=0; x < n; x++) {
                sum += ori[x] * Polytope_get_a(P, i, x);
            }
            distance[i] = Polytope_get_b(P, i) - sum;
#ifdef DEBUG
            //printf("%0.*f ", FLOATWIDTH, distance[i]);
#endif
        }
#ifdef DEBUG        
        //printf("\n");
#endif
     
        
        //check if ori in polytope
        for(i = 0; i < m; i++) {
            // sum = (A * ori)_i > b_i (constraint not satisfied)
            if(distance[i] < 0) {
                // tm[i] = row_i(A)*T*row_i(A)^t
#ifdef DEBUG
                //printf("%d in LOOP 1, distance[%d] = %0.15f < 0\n", i, i, distance[i]);
#endif
                tm[i] = 0;
                for (int j = 0; j < n; j++){
                    FT tmi_tmp = 0;
                    for (int k = 0; k < n; k++){
                        tmi_tmp += T[j*n + k] * Polytope_get_a(P, i, k);
                    }
                    tm[i] += Polytope_get_a(P, i, j) * tmi_tmp;
                }
                break;
            }
        }
        
        // check if small ellipsoid is contained in polytope
        if (i == m) {
            for (i = 0; i < m; i++){
#ifdef DEBUG
                //printf("%d in LOOP 2\n", i);
#endif
                tm[i] = 0;
                for (int j = 0; j < n; j++){
                    FT tmi_tmp = 0;
                    for (int k = 0; k < n; k++){
                        tmi_tmp += T[j*n + k] * Polytope_get_a(P, i, k);
                    }
                    tm[i] += Polytope_get_a(P, i, j) * tmi_tmp;
                }
  
                if (c3 * distance[i] * distance[i] - tm[i] < 0){
#ifdef DEBUG
                    //printf("%0.15f * %0.15f * %0.15f - %0.15f = %0.15f < 0\n", c3, distance[i], distance[i], tm[i], c3*distance[i]*distance[i]-tm[i]);
#endif
                    break;
                }
                else {
#ifdef DEBUG
                    //printf("%0.15f * %0.15f * %0.15f - %0.15f = %0.15f >= 0\n", c3, distance[i], distance[i], tm[i], c3*distance[i]*distance[i]-tm[i]);
#endif
                }
            }
        }

#ifdef PRINT_TMI
        //printf("ROUND %d\n", counter);
        for (int j = 0; j < m; j++){
            printf("%0.*f ", FLOATWIDTH, tm[j]);
        }
        printf("\n");
#endif

        //terminate if E satisfies the two criteria 
        if (i == m){
            break;
        }

        // else update ellipsoid (ori, T)
        for (int k = 0; k < n; k++){
            t[k] = 0;
            for (int j = 0; j < n; j++){
                t[k] += T[k*n+j] * Polytope_get_a(P, i, j);
            }
            if (tm[i] <= 0){

#ifdef DEBUG
                printf("tmi <= 0 for i = %d\nprinting tm\n", i);
                for (int l = 0; l < n; l++){
                    printf("%0.*f ", FLOATWIDTH, tm[i]);
                }
                printf("\n");
#endif
            }
            t[k] /= sqrt(tm[i]);
        }
        for (int k = 0; k < n; k++){
            ori[k] -= t[k] * c2;
#ifdef PRINT_T
            printf("%0.*f ", FLOATWIDTH, ori[k]);
#endif      
        }

#ifdef PRINT_T
        printf("\n\n");
#endif  
        
        for (int k = 0; k < n; k++){
            for (int j = 0; j < n; j++){
                T[k*n + j] = c1 * (T[k*n + j] - c4 * (t[k] * t[j]));
#ifdef PRINT_T
                printf("%0.*f ", FLOATWIDTH, T[k*n + j]);
#endif
            }
#ifdef PRINT_T
            printf("\n");
#endif
            
        }

#ifdef PRINT_T
        printf("\n");
#endif  

        
        
    }
    
    //apply affine transformation in-place on Poly
    FT *Trans = (FT *) calloc(n*n, sizeof(FT));
    int err = cholesky(T, Trans, n);
    if (err > 0){
        printf("The input polytope is degenerate or non-existant and the volume is 0.\n");
        exit(1);		
    }
    
    
    // b = beta_r * (b - A * ori);
    for (int i = 0; i < m; i++){
        Polytope_set_b(*Q, i, beta_r * distance[i]);
    }
    // A = A * Trans;
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            FT sum = 0;
            for (int k = 0; k < n; k++){
                sum += Polytope_get_a(P, i, k) * Trans[k*n + j];
            }
            Polytope_set_a(*Q, i, j, sum);
        }
    }

    /*
      don't know what that's supposed to do

    rowvec exp(n);
    exp.ones();
    for (int i = 0; i < n; i++){
        B[i] = b / A.col(i);
        Ai[i] = A / (A.col(i) * exp);
    }
    */

    *det = 1;
    for (int i = 0; i < n; i++){
        *det *= Trans[i*n+i];
    }
    *det /= pow(beta_r, n);

    free(T);
    free(distance);
    free(tm);
    free(t);
        
#ifdef DEBUG_MSG
    printf("Final ellipsoid\n");
    printf("\nT:\n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            printf("%0.*f ", FLOATWIDTH, T[i*n + j]);
        }
        printf("\n");
    }
    printf("\ncenter:\n");
    for (int i = 0; i < n; i++){
        printf("%0.*f ", FLOATWIDTH, ori[i]);
    }
    printf("\n");
    printf("\nTrans:\n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            printf("%0.*f ",FLOATWIDTH, Trans[i*n + j]);
        }
        printf("\n");
    }
    printf("\nTransformed Poly:\n");
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            printf("%0.*f ", FLOATWIDTH, Polytope_get_a(*Q, i, j));
        }
        printf("| %0.*f\n", FLOATWIDTH, Polytope_get_b(*Q, i));
    }

    printf("\n");

    printf("\nDeterminant:\n%0.*f\n", FLOATWIDTH, *det);
    
    printf("\nNumber of iterations of shallow beta-cut: %d\n", counter);

    printf("\n^^^^^^^^^^^^^^^^^ END HIGHVOLUMES\n");
    
#endif

    
}



//------------------------------
void preprocess_opcount(Polytope *P, FT R2, FT *ori, Polytope **Q, FT *det, int *iterations, int *loopone, int *looptwo, int *breakcond){

    int n = P->n;
    int m = P->m;

    *Q = Polytope_new(n, m);

    //------------------------------
    double beta_r = 2*n; 
    // (0,1,0,0)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //------------------------------
    double c1 = (2 * pow(n, 2) + pow(1 - n / beta_r, 2)) * (1 - 1.0 / pow(beta_r, 2)) / (2 * pow(n, 2) - 2);
    // (4,7,3,0)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //------------------------------
    double c2 = (1 - n / beta_r) / (n + 1);
    // (2, 0, 2, 0)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //------------------------------
    double c3 = beta_r * beta_r;
    // (0,1,0,0)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //------------------------------
    double c4 = 2 * c2 / (1 - 1.0 / beta_r);
    // (1, 1, 2, 0)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    // initialize T to diag(R2)
    // T is out initial guess to the ellipsoid around Poly
    FT *T = (FT *) calloc(n*n, sizeof(FT));
    for (int i = 0; i < n; i++){
        T[i * n + i] = R2;
    }

    FT *distance = (FT *) calloc(m, sizeof(FT));
    FT *tm = (FT *) calloc(m, sizeof(FT));
    FT *t = (FT *) malloc(n*sizeof(FT));

    // update T until x^t*T*x <= 1

    *iterations = 0;
    *loopone = 0;
    *looptwo = 0;
    *breakcond = 0;
    //------------------------------
    while (++(*iterations) > 0){
        int i;

        // distance = b - A * ori
        //------------------------------
        for (i = 0; i < m; i++){
            (*loopone)++;

            //------------------------------
            FT *Ai = Polytope_get_Ai(P, i);
            distance[i] = Polytope_get_b(P, i);
            distance[i] -= dotProduct(Ai, ori, n);
            // (n+1, n)
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            if (distance[i] < 0){
                // tm[i] = row_i(A)*T*row_i(A)^t
                (*breakcond)++;

                //------------------------------
                tm[i] = 0;
                for (int j = 0; j < n; j++){
                    FT tmi_tmp = dotProduct(Ai, &T[j*n], n);
                    tm[i] += Polytope_get_a(P, i, j) * tmi_tmp;
                }
                // (n^2 + n, n^2 + n)
                //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                break;
            }
        }
        // loopone * (n+1, n) + breakcond *(n^2 + n, n^2 + n)
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        // check if small ellipsoid is contained in polytope
        if (i == m) {
            //------------------------------
            for (i = 0; i < m; i++){
                (*looptwo)++;
                //------------------------------
                tm[i] = 0;
                for (int j = 0; j < n; j++){
                    FT tmi_tmp = dotProduct(Polytope_get_Ai(P, i), &T[j*n], n);
                    tm[i] += Polytope_get_a(P, i, j) * tmi_tmp;
                }
                // (n^2 + n, n^2 + n)
                //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                //------------------------------
                if (c3 * distance[i] * distance[i] - tm[i] < 0){
                    break;
                }
                // (1, 2)
                //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            }
            // looptwo * (n^2 + n + 1, n^2 + n + 2)
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        }


        //terminate if E satisfies the two criteria 
        if (i == m){
            break;
        }

        // else update ellipsoid (ori, T)
        //------------------------------
        for (int k = 0; k < n; k++){
            t[k] = dotProduct(Polytope_get_Ai(P, i), &T[k*n], n);
            t[k] /= sqrt(tm[i]);
        }
        // n * (n, n, 1, 1)
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        //------------------------------
        for (int k = 0; k < n; k++){
            ori[k] -= t[k] * c2;   
        }
        // n * (1, 1)
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        //------------------------------
        for (int k = 0; k < n; k++){
            for (int j = 0; j < n; j++){
                T[k*n + j] = c1 * (T[k*n + j] - c4 * (t[k] * t[j]));
            }   
        }
        // n^2 * (1, 3)
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    }
    // iterations * (2n^2 + n, 4n^2 + n, n, n) +
    // loopone * (n+1, n) + 
    // breakcond *(n^2 + n, n^2 + n) +
    // looptwo * (n^2 + n + 1, n^2 + n + 2)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    //apply affine transformation in-place on Poly
    FT *Trans = (FT *) calloc(n*n, sizeof(FT));
    //------------------------------
    int err = cholesky(T, Trans, n);
    // ((n^4 - n^3 + 2n^2 + 4n)/6,
    //  (n^4 - 4n^3 + 8n^2 - 5n)/6,
    //  (n^2 - n)/2,
    //  n)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (err > 0){
        printf("The input polytope is degenerate or non-existant and the volume is 0.\n");
        exit(1);		
    }
    
    
    // b = beta_r * (b - A * ori);
    //------------------------------
    for (int i = 0; i < m; i++){
        Polytope_set_b(*Q, i, beta_r * distance[i]);
    }
    // (0, m, 0, 0)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    // A = A * Trans;
    //------------------------------
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            FT sum = 0;
            for (int k = j; k < n; k++){
                sum += Polytope_get_a(P, i, k) * Trans[k*n + j];
            }
            Polytope_set_a(*Q, i, j, sum);
        }
    }
    // (mn(n+1)/2, mn(n+1)/2, 0, 0)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    //------------------------------
    *det = 1;
    for (int i = 0; i < n; i++){
        *det *= Trans[i*n+i];
    }
    *det /= pow(beta_r, n);
    // (0, 2n-1, 1, 0)
    // question: pow(x, n) is n-1 mults?
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    
    free(T);
    free(distance);
    free(tm);
    free(t);
}
// adds: iterations * (2n^2 + n) +
//       loopone * (n+1) +
//       breakcond * (n^2 + n) + 
//       looptwo * (n^2+n+1) +
//       (n^3 + 6n^2 + 5n)/6 + mn(n+1)/2 + 7
// mults: iterations * (4n^2 + n) +
//        loopone * n +
//        breakcond * (n^2 + n) + 
//        looptwo * (n^2+n+2) +
//        (n^3 + 3n^2 - 4n)/6 + mn(n+1)/2 + 2n + m + 9
// divs: iterations * n +
//       (n^2 + n)/2 + 8
// sqrts: iterations * n +
//        n
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    // = ((n^3 + 6n^2 + 5n)/6,
    //    (n^3 + 3n^2 - 4n)/6,
    //    (n^2 + n)/2,
    //    n)
