
#include <stdlib.h>
#include <string.h>
#include <glpk.h>
//#include "beta_cut.h"
#include "preprocess.h"



void init_ellipsoid(const Polytope *Pol, FT *R2, FT **Ori){
    int n = Pol->n;
    int m = Pol->m;
    
    *R2 = 0;
    *Ori = (FT *) calloc(n, sizeof(FT));

    //init GLPK
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, m);
    glp_add_cols(lp, n);

    //disable msg output
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR;

    //load constraints
    // MB: is this 1-based indexing due to glp?
    int *ind = (int *) malloc ((n+1) * sizeof(int));
    FT *val = (FT *) malloc((n+1) * sizeof(FT));
    for (int i = 1; i < m + 1; i++){
        for (int j = 1; j < n + 1; j++){
            ind[j] = j;
            val[j] = Polytope_get_a(Pol, i-1, j-1);
        }
        
        glp_set_mat_row(lp, i, n, ind, val);
        // note: A[i][n] = b[i]
        glp_set_row_bnds(lp, i, GLP_UP, 0, Polytope_get_b(Pol, i-1));
        
    }
    free(ind);
    free(val);
    
    for (int i = 1; i < n + 1; i++){
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
    }

    
    //get max/min bounds of polytope in each of the n dimensions
    for (int i = 0; i < n; i++){
        FT max, min;
        for (int j = 0; j < n; j++) {
            glp_set_obj_coef(lp, j + 1, 0);
        }

        // max bound
        glp_set_obj_coef(lp, i + 1, 1);
        glp_simplex(lp, &parm);
        max = glp_get_obj_val(lp);
        for (int j = 0; j < n; j++) {
            (*Ori)[j] += glp_get_col_prim(lp, j + 1);
        }

        // min bound
        glp_set_obj_coef(lp, i + 1, -1);
        glp_simplex(lp, &parm);
        min = -glp_get_obj_val(lp);
        for (int j = 0; j < n; j++)
            (*Ori)[j] += glp_get_col_prim(lp, j + 1);

        // note: sqrt(R2) will be diameter of bounding box 
        *R2 += (max - min) * (max - min);
    }
    for (int i = 0; i < n; i++){
        (*Ori)[i] /= (2 * n);
    }
	
    glp_delete_prob(lp);
}



// TODO: put in own function
void Polytope_bounding_ref(const void *B, FT *R2, FT **Ori){
    const Polytope *P = (Polytope *) B;
    init_ellipsoid(P, R2, Ori);
}


void Ellipsoid_bounding_ref(const void *B, FT *R2, FT **ori){

    const Ellipsoid *E = (Ellipsoid *) B;
    int n = E->n;
    
    // center is trivially given
    *ori = E->a;
    
    // maximize each of the 2n linear functions e_i, -e_i
    // we could also compute eigenvalues of E but this seems more difficult
    // note we can ignore center a for this
    // note the linear function e_i is maximized by sqrt(T_ii)
    // note we choose R2 = sum_{i} (max e_i - min e_i)^2 as for polytopes
    Matrix *T = Matrix_new(n,n);

    // once ellipsoid has matrix members this copying won't be necessary
    // but for now we do it, as we want to be 32byte aligned in as many functions as possible
    Matrix *A = Matrix_new(n,n);
    for (int i = 0; i < n; i++){
        FT *Ai = Matrix_get_row(A, i);
        FT *Ei = Ellipsoid_get_Ai(E, i);
        for (int j = 0; j < n; j++){
            Ai[j] = Ei[j];
        }
    }

    
    Matrix_invert_pdsym(A, T);
    for (int i = 0; i < n; i++){
        *R2 += 4*Matrix_get_row(T, i)[i];
    }  
}
