#include "test_helpers.hpp"


/**
 * \brief compare the two polytopes elementwise and return the 2-frobenius norm of the difference matrix of A and the 2-norm of the difference of b
 * we consider EPS as 0
 **/
std::pair<FT, FT> matrix_diff(Polytope *P, vol::Polyvest_p *Q){

    int n = P->n;
    int m = P->m;

    std::pair<FT, FT> res = make_pair(0.0, 0.0);
    
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            FT a1 = Polytope_get_a(P, i, j);
            FT a2 = Q->A(i, j);
            FT a_diff = a1 - a2;
            a_diff *= a_diff;
            if (a_diff > EPS){
                res.first += a_diff;
            }
            if (!std::isfinite(a1) || !std::isfinite(a2)){
                res = std::make_pair(-1, -1);
                break;
            }
        }
        FT b1 = Polytope_get_b(P, i);
        FT b2 = Q->b(i);
        FT b_diff = b1 - b2; 
        b_diff *= b_diff;
        if (b_diff > EPS){
            res.second += b_diff;
        }        
        if (!std::isfinite(b1) || !std::isfinite(b2)){
            res = std::make_pair(-1, -1);
            break;
        }
    }

    return res;
    
}




FT frobenius(FT *A, FT *B, int d1, int d2) {
    //Polytope *P, vol::Polyvest_p *Q){

    FT res = 0.0;
    
    for (int i = 0; i < d1; i++){
        for (int j = 0; j < d2; j++){
            FT a1 = A[i*d2+j];
            FT a2 = B[i*d2+j];
            FT a_diff = a1 - a2;
            a_diff *= a_diff;
            if (a_diff > EPS){
                res += a_diff;
            }
            if (!std::isfinite(a1) || !std::isfinite(a2)){
                res = -1;
                break;
            }
        }
    }

    return res;
    
}




bool polytope_contains_scaled_ball(Polytope *P){

    // setup scaled unit ball for ellipsoid_inside_poly method
    int n = P->n;

    FT *c = (FT *) calloc(n, sizeof(FT));
    FT *E = (FT *) calloc(n*n, sizeof(FT));
    for (int i = 0; i < n; i++){
        // note beta = 2n^{-1}
        E[i*n+i] = 1.0/(2*n);
    }

    return ellipsoid_inside_poly(P, E, c);
    
}


bool ellipsoid_inside_poly(Polytope *P, FT *E, FT *c){

    // we maximize each linear constraint of P inside the ellipsoid and see if it is satisfied
    // note that the linear function f(x) = a.transpose() * x is maximized in (E,c) by a.transpose() * c + sqrt(a.transpose() * E * a)
    // note this wasn't obvious to me but i read it in the shallow beta-cut paper P. 69

    int n = P->n;
    int m = P->m;
    
    for (int i = 0; i < m; i++){
        
        // bi <- b[i] - Ai * c
        FT bi = Polytope_get_b(P, i);
        for (int j = 0; j < n; j++){
            bi -= Polytope_get_a(P, i, j) * c[j];
        }

        // at_e_a <- Ai * E * Ai.transpose()
        FT at_e_a = 0;
        for (int j = 0; j < n; j++){
            FT e_a = 0;
            for (int k = 0; k < n; k++){
                e_a += E[j*n+k] * Polytope_get_a(P, i, k);
            }
            at_e_a += Polytope_get_a(P, i, j) * e_a;
        }

        //  Ai * c + sqrt(Ai * E * Ai.transpose()) > bi 
        if (at_e_a > bi*bi){
            return false;
        }
        
    }
    
    return true;
    
}



void polyvest_convert(Polytope *P, vol::Polyvest_p *Q){

    int n = P->n;
    int m = P->m;  

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            Q->matA(Polytope_get_a(P, i, j), i, j);
        }
        Q->vecb(Polytope_get_b(P, i), i);
    }
  
}


int read_polyvest_p(string filename, Polytope **P){

    ifstream file;
    file.open(filename);


    if (!file.is_open()){
        printf("failed to read polytope");
        return 1;
    }

    int n, m;
    file >> m >> n;

    *P = Polytope_new(n, m);

    FT num;
    for (int i = 0; i < m; i++){
        file >> num;
        Polytope_set_b(*P, i, num);
        for (int j = 0; j < n; j++){
            file >> num;
            Polytope_set_a(*P, i, j, num);
        }
    }

    return 0;
}
