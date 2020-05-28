#include "polytope.h"
#include <glpk.h>
#include <omp.h>

glp_prob* get_lp(const int n, const int m) {
    //init GLPK
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, m);
    glp_add_cols(lp, n);

    return lp; 
}

typedef glp_prob* (*X_get_lp_f_t) (const void*);

glp_prob* Polytope_get_lp(const void* o) {
    Polytope* P = (Polytope*)o;
    glp_prob* lp = get_lp(P->n,P->m);
    const int n = P->n;
    const int m = P->m;

    //load constraints
    // MB: is this 1-based indexing due to glp?
    int *ind = (int *) malloc ((n+1) * sizeof(int));
    FT *val = (FT *) malloc((n+1) * sizeof(FT));
    for (int i = 1; i < m + 1; i++){
        for (int j = 1; j < n + 1; j++){
            ind[j] = j;
            val[j] = Polytope_get_a(P, i-1, j-1);
        }
        
        glp_set_mat_row(lp, i, n, ind, val);
        // note: A[i][n] = b[i]
        glp_set_row_bnds(lp, i, GLP_UP, 0, Polytope_get_b(P, i-1));
        
    }
    free(ind);
    free(val);

    // constraints on variables: free
    for (int i = 1; i < n + 1; i++){
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
    }

    return lp;
}

glp_prob* PolytopeT_get_lp(const void* o) {
    PolytopeT* P = (PolytopeT*)o;
    glp_prob* lp = get_lp(P->n,P->m);
    const int n = P->n;
    const int m = P->m;

    //load constraints
    // MB: is this 1-based indexing due to glp?
    int *ind = (int *) malloc ((n+1) * sizeof(int));
    FT *val = (FT *) malloc((n+1) * sizeof(FT));
    for (int i = 1; i < m + 1; i++){
        for (int j = 1; j < n + 1; j++){
            ind[j] = j;
            val[j] = PolytopeT_get_a(P, i-1, j-1);
        }
        
        glp_set_mat_row(lp, i, n, ind, val);
        // note: A[i][n] = b[i]
        glp_set_row_bnds(lp, i, GLP_UP, 0, PolytopeT_get_b(P, i-1));
        
    }
    free(ind);
    free(val);

    // constraints on variables: free
    for (int i = 1; i < n + 1; i++){
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
    }

    return lp;
}


glp_prob* PolytopeCSC_get_lp(const void* o) {
    const PolytopeCSC* P = (const PolytopeCSC*)o;
    glp_prob* lp = get_lp(P->n,P->m);
    const int n = P->n;
    const int m = P->m;

    //load constraints
    // MB: is this 1-based indexing due to glp?
    int *ind = (int *) malloc ((n+1) * sizeof(int));
    for (int j = 1; j < n + 1; j++){
        ind[j] = j;
    }

    FT *val = (FT *) malloc((n+1) * sizeof(FT));
    for (int i = 1; i < m + 1; i++){
        PolytopeCSC_get_Ai(P, i-1, val+1);        
        glp_set_mat_row(lp, i, n, ind, val);
        // note: A[i][n] = b[i]
        glp_set_row_bnds(lp, i, GLP_UP, 0, P->b[i-1]);
        
    }
    free(ind);
    free(val);

    // constraints on variables: free
    for (int i = 1; i < n + 1; i++){
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
    }

    return lp;
}



void LP_to_boundingSphere(const int n, X_get_lp_f_t x_get_lp, const void* o, FT *R2, FT *Ori) {
    *R2 = 0;
    for (int j = 0; j < n; j++) {Ori[j] = 0;}
 
    #pragma omp parallel
    {
       int tid = omp_get_thread_num();
       if(volumeVerbose>=2) {printf("Hello from %d\n",tid); }
       
       FT R2_local = 0;
       FT Ori_local[n];
       for(int i=0;i<n;i++) {Ori_local[i]=0;}
       
       // generate lp:
       glp_prob* lp = x_get_lp(o);

       //disable msg output
       glp_smcp parm;
       glp_init_smcp(&parm);
       parm.msg_lev = GLP_MSG_ERR;

       //get max/min bounds in each of the n dimensions
       #pragma omp for schedule(dynamic)
       for (int i = 0; i < n; i++){
	   if(volumeVerbose>=2) {printf(" - thread %d does simplex for dimension %d.\n",tid,i); }
           FT max, min;
           for (int j = 0; j < n; j++) {
               glp_set_obj_coef(lp, j + 1, 0);
           }

           // max bound
           glp_set_obj_coef(lp, i + 1, 1);
           glp_simplex(lp, &parm);
           max = glp_get_obj_val(lp);
           for (int j = 0; j < n; j++) {
               Ori[j] += glp_get_col_prim(lp, j + 1);
           }

           // min bound
           glp_set_obj_coef(lp, i + 1, -1);
           glp_simplex(lp, &parm);
           min = -glp_get_obj_val(lp);
           for (int j = 0; j < n; j++)
               Ori_local[j] += glp_get_col_prim(lp, j + 1);

           // note: sqrt(R2) will be diameter of bounding box 
           R2_local += (max - min) * (max - min);
       }
       glp_delete_prob(lp);
       
       // write back:
       #pragma omp critical
       {
          *R2 += R2_local;
	  for(int i=0;i<n;i++) {Ori[i] += Ori_local[i];}
       } 
    }
   
    for (int i = 0; i < n; i++){
        Ori[i] /= (2 * n);
    }
	
}

void Polytope_bounding_ref(const void *B, FT *R2, FT *Ori){
    const Polytope *P = (Polytope *) B;
    LP_to_boundingSphere(P->n, &Polytope_get_lp, P, R2, Ori);
}

void PolytopeT_bounding_ref(const void *B, FT *R2, FT *Ori){
    const PolytopeT *P = (PolytopeT *) B;
    LP_to_boundingSphere(P->n, &PolytopeT_get_lp, P, R2, Ori);
}

void PolytopeCSC_bounding_ref(const void *o, FT *r, FT *ori){
    const PolytopeCSC *P = (PolytopeCSC *) o;
    LP_to_boundingSphere(P->n, &PolytopeCSC_get_lp, P, r, ori);
}

