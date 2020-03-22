#include "beta_cut.h"



void Preprocess(double beta_r, Polytope *P){

    // MB: maybe implement this function as in PolyVest
    //checkHPs();

    int n = P->n;
    int m = P->m;

    double c3 = beta_r * beta_r;
    double c1 = (2 * n*n + (1-n/beta_r)*(1-n/beta_r)) * (1 - 1.0 / c3) / (2 * n*n - 2);
    double c2 = (1 - n / beta_r) / (n + 1);
    double c4 = 2 * c2 / (1 - 1.0 / beta_r);

    //init E(R2I, 0), T = R2I, ori = 0.
    
    FT R2;
    FT *ori;
    initEllipsoid(P, &R2, &ori);

    // initialize T to diag(R2)
    FT *T = (FT *) calloc(n*n, sizeof(FT));
    for (int i = 0; i < n; i++){
        T[i * n + i] = R2;
    }

    FT *distance = (FT *) calloc(m, sizeof(FT));
    FT *tm = (FT *) calloc(m, sizeof(FT));

    int counter = 0;
    while (++counter > 0){
        int i;
		
        //check if ori in polytope
        for(i = 0; i < m; i++) {
            FT sum = 0;
            for(int x=0; x < n; x++) {
                sum += ori[x] * Polytope_get_a(P, i, x);
            }
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

    
    //apply affine transformation
    //mat Trans = chol(T);

    // TODO from here!
    FT *Trans = (FT *) malloc(n*n*sizeof(FT));
    int err = cholesky(&Trans, T);
    if (err > 0){
        printf("The input polytope is degenerated or non-existed and the volume is 0.\n");
        exit(1);		
    }

    /*cout << Trans << endl;*/
    b = beta_r * (b - A * ori);
    A = A * Trans.t();

    if (!msg_off) cout << "The number of iterations: " << counter << endl;

    rowvec exp(n);
    exp.ones();
    for (int i = 0; i < n; i++){
        B[i] = b / A.col(i);
        Ai[i] = A / (A.col(i) * exp);
    }
	
    determinant = det(Trans) / pow(beta_r, n);
}
