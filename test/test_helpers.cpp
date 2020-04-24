#include "test_helpers.hpp"
#include <glpk.h>

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
            FT a_diff = abs(a1 - a2);
            //a_diff *= a_diff;
            //cout << a_diff << std::endl;
            if (a_diff > EPS){
                res.first += a_diff;
            }
            if (!std::isfinite(a1) || !std::isfinite(a2)){
                res = std::make_pair(-1, -1);
                return res;
            }
        }
        FT b1 = Polytope_get_b(P, i);
        FT b2 = Q->b(i);
        FT b_diff = abs(b1 - b2); 
        //b_diff *= b_diff;
        if (b_diff > EPS){
            res.second += b_diff;
        }
        if (!std::isfinite(b1) || !std::isfinite(b2)){
            res = std::make_pair(-1, -1);
            return res;
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
                return res;
            }
        }
    }

    return res;
    
}




bool polytope_in_unit_ball(Polytope *P){


    int n = P->n;
    int m = P->m;

    //init GLPK
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, m);
    glp_add_cols(lp, n);

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR;

    int *ind = new int[n + 1];
    for (int j = 1; j < n + 1; j++){
        // setup for constraint loading
        ind[j] = j;
        // all variables are free
        glp_set_col_bnds(lp, j, GLP_FR, 0, 0);
        // only check feasibility
        glp_set_obj_coef(lp, j, 1);
    }
    // load polytope constraints
    for (int i = 1; i < m+1; i++){
        glp_set_mat_row(lp, i, n, ind, Polytope_get_Ai(P, i-1) - 1);
    }
    
    
    // use this bit-magic to compute all subsets of m of size n
    // this avoids computing all subsets of m and filtering out the relevant ones
    unsigned int v = (1 << n) - 1; // start with lexicographically smallest permutation of n ones
    unsigned int w = 0;

    while (v < (1 << m)){

        
        // check the vertex specified by v
        for (int i = 0; i < m; i++){
            if (v & (1 << i)){ // force constraint 
                glp_set_row_bnds(lp, i+1, GLP_FX, Polytope_get_b(P, i), Polytope_get_b(P, i));
            }
            else if (w & (1 << i)){ // constraint was only forced in last iteration
                glp_set_row_bnds(lp, i+1, GLP_UP, 0, Polytope_get_b(P, i));
            }
        }
        w = v;

        glp_simplex(lp, &parm);

        // if problem is infeasible the n constraints are not linearly independent
        // if n constraints are dependent and the problem is feasible, the following check still needs to hold as we are inside the polytope, although not on a vertex
        if (glp_get_status(lp) == GLP_OPT || glp_get_status(lp) == GLP_FEAS){
            // check that point in unit ball
            double sum = 0;
            //std::cout << "vertex: ( ";
            for (int i = 1; i < n+1; i++){
                double val = glp_get_col_prim(lp, i);
                //std::cout << val << " ";
                sum += val*val;
            }
            //std::cout << ")\n";
            if (sum > 1){
                return false;
            }
        }
                
        // compute lexicographically next permutation
        unsigned int t = v | (v - 1); // t gets v's least significant 0 bits set to 1
        // Next set to 1 the most significant bit to change, 
        // set to 0 the least significant ones, and add the necessary 1 bits.
        v = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
    }

    glp_delete_prob(lp);
    delete []ind;
    
    return true;
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
        
        // b[i] - Ai * c
        FT bi_sub_Ai_dot_c = Polytope_get_b(P, i);
        for (int j = 0; j < n; j++){
            bi_sub_Ai_dot_c -= Polytope_get_a(P, i, j) * c[j];
        }

        // at_e_a <- Ai.transpose() * E * Ai
        FT at_e_a = 0;
        for (int j = 0; j < n; j++){
            FT e_a = 0;
            for (int k = 0; k < n; k++){
                e_a += E[j*n+k] * Polytope_get_a(P, i, k);
            }
            at_e_a += Polytope_get_a(P, i, j) * e_a;
        }

        // sqrt(Ai.transpose() * E * Ai) > bi - Ai * c
        if (at_e_a > bi_sub_Ai_dot_c * bi_sub_Ai_dot_c){
            return false;
        }
        
    }
    
    return true;
    
}
