#include "polytopeCSC.h"


Body_T PolytopeCSC_T = {
	.print = PolytopeCSC_print,
        .free = PolytopeCSC_free,
        .clone = PolytopeCSC_clone,
        .inside = PolytopeCSC_inside_ref,
        .intersect = PolytopeCSC_intersect_ref,
        .intersectCoord = PolytopeCSC_intersectCoord_ref,
	.cacheAlloc = PolytopeCSC_cacheAlloc_ref,
	.cacheReset = PolytopeCSC_cacheReset_ref,
	.cacheUpdateCoord = PolytopeCSC_cacheUpdateCoord_ref,
        .shallowCutOracle = PolytopeCSC_shallowCutOracle_ref,
	.transform = PolytopeCSC_transform_ref,
        .boundingSphere = PolytopeCSC_bounding_ref
};





void PolytopeCSC_print(const void* o) {
   const PolytopeCSC *p = (PolytopeCSC *)o;
   printf("PolytopeCSC: n=%d, m=%d\n",p->n,p->m);
   printf("Printing raw data!\n");

   printf("A     row\n");
   for (int i = 0; i < p->n; i++){
       printf("-------- Col %d\n", i);
       for (int j = p->col_start[i]; j < p->col_start[i+1]; j++){
           printf("%.3f %d\n", p->A[j], p->row_idx[j]);
       }
   }

   
   printf("\n\nPrinting transpose!\n");
   
   for(int j = 0; j < p->n; j++) {
       int idx = 0;
       for (int jj = p->col_start[j]; jj < p->col_start[j+1], p->row_idx[jj] > -1; jj++){
           while (idx < p->row_idx[jj]){
               printf(" %.3f", 0.);
               idx++;
           }
           printf(" %.3f", p->A[jj]);
           idx++;
       }
       for (; idx < p->m; idx++){
           printf(" %.3f", 0.);
       }
       printf("\n");
   }

   for (int i = 0; i < p->m; i++){
       printf("%.3f\n", p->b[i]);
   }
}


void *PolytopeCSC_clone(const void *o){
    PolytopeCSC *p = (PolytopeCSC *) o;
    PolytopeCSC *p_clone = (PolytopeCSC *) malloc(sizeof(PolytopeCSC));
    
    p_clone->n = p->n;
    p_clone->m = p->m;

    p_clone->b = (FT *) aligned_alloc(32, p_clone->m * sizeof(FT));
    memcpy(p_clone->b, p->b, p_clone->m * sizeof(FT));
    
    p_clone->col_start = (int *) aligned_alloc(32, (p_clone->n+1)*sizeof(int));
    memcpy(p_clone->col_start, p->col_start, (p_clone->n+1)*sizeof(int));

    int ndata = p_clone->col_start[p_clone->n];
    
    p_clone->A = (FT *) aligned_alloc(32, ndata * sizeof(FT));
    memcpy(p_clone->A, p->A, ndata * sizeof(FT));

    p_clone->row_idx = (int *) aligned_alloc(32, ndata * sizeof(int));
    memcpy(p_clone->row_idx, p->row_idx, ndata * sizeof(int));

    return (void *) p_clone;
}

void PolytopeCSC_free(const void *o){
    PolytopeCSC *p = (PolytopeCSC *) o;
    free(p->A);
    free(p->b);
    free(p->row_idx);
    free(p->col_start);
    free(p);
}


bool PolytopeCSC_inside_ref(const void *o, const FT *x){
    PolytopeCSC *p = (PolytopeCSC *) o;

    int n = p->n;
    int m = p->m;
    
    // check all inequalities
    // m accumulators, one per constraint
    FT *acc = (FT *) calloc(m, sizeof(FT));

    // walk through p and update at each valid entry the corresponding acc
    for (int i = 0; i < n; i++){
        FT xi = x[i];
        for (int j = p->col_start[i]; j < p->col_start[i+1], p->row_idx[j] > -1; j++){
            acc[p->row_idx[j]] += xi * p->A[j];
        }
    }

    // compare accs with all bs
    bool res = true;
    for (int i = 0; i < m; i++){
        if (acc[i] > p->b[i]){
            res = false;
        }
    }

    free(acc);
    return res;
    
}


void PolytopeCSC_intersect_ref(const void *o, const FT *x, const FT *dir, FT *t0, FT *t1){

    PolytopeCSC *p = (PolytopeCSC *) o;

    // make vector aligned in case we want to vectorize this at some point
    FT *dotd = dotproduct_store_d;
    FT *dotx = dotproduct_store_x;
    memset(dotd, 0, p->m*sizeof(FT));
    memset(dotx, 0, p->m*sizeof(FT));

    for (int i = 0; i < p->n; i++){
        FT xi = x[i];
        FT di = dir[i];
        for (int j = p->col_start[i]; j < p->col_start[i+1], p->row_idx[j] > -1; j++){
            dotx[p->row_idx[j]] += p->A[j] * xi;
            dotd[p->row_idx[j]] += p->A[j] * di;
        }
    }

    *t0 = -FT_MAX;
    *t1 = FT_MAX;
    
    FT *b = p->b;
    for (int i = 0; i < p->m; i++){
        // only consider constraints that are non-collinear with dir
        if (dotd[i] > FT_EPS || -dotd[i] > FT_EPS) {
            FT t = (b[i] - dotx[i]) / dotd[i];

            if (dotd[i] < 0.0) {
                *t0 = (*t0 > t) ? *t0 : t;
            }
            else {
                *t1 = (*t1 < t) ? *t1 : t;
            }
            
        }
    }

}


void PolytopeCSC_mvm(const PolytopeCSC *p, const FT *x, FT *res){

    int m = p->m;
    int n = p->n;
    memset(res, 0, m*sizeof(FT));

    for (int i = 0; i < n; i++){
        FT xi = x[i];
        for (int j = p->col_start[i]; j < p->col_start[i+1], p->row_idx[j] > -1; j++){
            res[p->row_idx[j]] += p->A[j] * xi;
        }
    }    
}


void PolytopeCSC_intersectCoord_ref(const void *o, const FT *x, const int d, FT *t0, FT *t1, void *cache) {

    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    FT *Aix = (FT *) cache;
    FT *b = p->b;
    
    *t0 = -FT_MAX;
    *t1 = FT_MAX;

    // compute full MVM Ax here
    // note that this is as cheap as computing only one dot product Ai x!!
    FT *dotx = dotproduct_store_x;
    PolytopeCSC_mvm(p, x, dotx);
    
    for (int i = p->col_start[d]; i < p->col_start[d+1], p->row_idx[i] > -1; i++){
        FT bi = b[p->row_idx[i]];
        FT dai = p->A[i];

        assert(dotx[p->row_idx[i]] == Aix[p->row_idx[i]] && "Cache must be accurate!");

        FT t = (bi - dotx[p->row_idx[i]]) / dai;

        if (dai < 0.0){
            *t0 = (*t0 > t) ? *t0 : t;
        }
        else {
            *t1 = (*t1 < t) ? *t1 : t;
        }
    }

}



int PolytopeCSC_cacheAlloc_ref(const void *o){
    // need to cache m dot products, one per inequality, i.e. m FT
    const PolytopeCSC *p = (PolytopeCSC *) o;
    return p->m * sizeof(FT);
}

void PolytopeCSC_cacheReset_ref(const void *o, const FT *x, void *cache){

    // compute dot product
    PolytopeCSC *p = (PolytopeCSC *) o;
    FT *c = (FT *) cache;
    memset(c, 0, p->m*sizeof(FT));

    for (int i = 0; i < p->n; i++){
        FT xi = x[i];
        for (int j = p->col_start[i]; j < p->col_start[i+1], p->row_idx[j] > -1; j++){
            c[p->row_idx[j]] += p->A[j] * xi;
        }
    }
}

void PolytopeCSC_cacheUpdateCoord_ref(const void *o, const int d, const FT dx, void *cache) {

    const PolytopeCSC *p = (PolytopeCSC *) o;
    FT *c = (FT *) cache;

    // only update with column d of A
    for (int i = p->col_start[d]; i < p->col_start[d+1], p->row_idx[i] > -1; i++){
        c[p->row_idx[i]] += dx * p->A[i];
    }
}




FT *PolytopeCSC_get_Ai(const PolytopeCSC *p, int row, FT *res){

    int n = p->n;
    
    // read row i into res
    for (int i = 0; i < n; i++){
            
        int j = p->col_start[i];
        while (j < p->col_start[i+1] &&
               p->row_idx[j] > -1 &&
               p->row_idx[j] < row){
            j++;
        }
        if (p->row_idx[j] == row){
            res[i] = p->A[j];
        }
        else {
            res[i] = 0.0;
        }
    } 
}


bool PolytopeCSC_shallowCutOracle_ref(const void* o, const Ellipsoid* E, FT* v, FT* c) {

    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;

    FT Ax[m];
    
    // check if center of ellipsoid in polytope
    FT *Ai = (FT *) aligned_alloc(32, n*sizeof(FT));
    for (int i = 0; i < m; i++){
        PolytopeCSC_get_Ai(p, i, Ai);
        FT bi = p->b[i];
        FT *x = E->a;
        Ax[i] = dotProduct(Ai, x, n);
        if (Ax[i] > bi){ // found one -> return (Ai, bi)
            for(int j=0;j<n;j++) {v[j] = Ai[j];}
            *c = bi;
            return true;
        }        
    }
    
    // compute full MVM A * E->a as this is as expensive as computing a single dot product
    FT *AEa = (FT *) aligned_alloc(32, m * sizeof(FT));
    PolytopeCSC_mvm(p, E->a, AEa);

    
    // this is costly and tedious
    // maybe we should do preprocessing on transformed body Polytope or PolytopeT and then transform back?
    const FT twon2 = 4.0*n*n;
    for (int i = 0; i < m; i++){
        PolytopeCSC_get_Ai(p, i, Ai);
        FT bi = p->b[i];
        FT AitTAi = 0; // could be useful to cache...
        for(int j=0;j<n;j++) {
            FT* Tj = Ellipsoid_get_Ti(E,j);
            FT TjAi = dotProduct(Tj,Ai,n);
            AitTAi += Ai[j] * TjAi;
        }
      
        FT diff = bi - Ax[i];
        if(AitTAi > diff*diff*twon2) { // found one -> return (Ai, bi)
            for(int j=0;j<n;j++) {v[j] = Ai[j];}
            *c = bi;
            return true;
        }
    }

    free(Ai);

    return false;

}



void PolytopeCSC_transform_ref(const void *o_in, void *o_out, const Matrix *L, const FT *a, const FT beta){

    const PolytopeCSC *p = (PolytopeCSC *) o_in;
    PolytopeCSC *p_out = (PolytopeCSC *) o_out;

    int n = p->n;
    int m = p->m;

    p_out->n = n;
    p_out->m = m;

    
    // col-major matrix to hold intermediate results of computations
    Matrix *M = Matrix_new(m, n);

    // compute A*L
    // go throught A in col major order
    for (int i = 0; i < n; i++){
        for (int j = p->col_start[i]; j < p->col_start[i+1], p->row_idx[j] > -1; j++){
            for (int k = 0; k <= i; k++){
                FT val = Matrix_get(M, k, p->row_idx[j]) + p->A[j] * Matrix_get(L, i, k);
                Matrix_set(M, k, p->row_idx[j], val);
            }
        }
    }

    // insert computed values into p_out
    // maybe not necessary to align this?
    p_out->col_start = (int *) aligned_alloc(32, (n+1)*sizeof(int));
    memset(p_out->col_start, 0, n+1);
    
    // first count nr of elements per column
    for (int i = 0; i < p_out->n; i++){
        p_out->col_start[i+1] = p_out->col_start[i];
        for (int j = 0; j < p_out->m; j++){
            if (Matrix_get(M, i, j) != 0){
                p_out->col_start[i+1]++;
            }
        }
        int cs = ceil_cache(p_out->col_start[i+1], sizeof(FT)) / 8;
        p_out->col_start[i+1] = cs;
    }

    
    // second fill values into A and row indices into row_idx
    p_out->A = (FT *) aligned_alloc(32, p_out->col_start[p_out->n] * sizeof(FT));
    p_out->row_idx = (int *) aligned_alloc(32, p_out->col_start[p_out->n] * sizeof(int));

    int idx = 0;
    for (int i = 0; i < p_out->n; i++){
        for (int j = 0; j < p_out->m; j++){
            if (Matrix_get(M, i, j) != 0){
                p_out->A[idx] = Matrix_get(M, i, j);
                p_out->row_idx[idx] = j;
                idx++;
            }
        }
        for (; idx < p_out->col_start[i+1]; idx++){
            p_out->row_idx[idx] = -1;
        }
    }

    // compute updated b and load into p_out
    p_out->b = (FT *) aligned_alloc(32, p_out->m * sizeof(FT));    
    PolytopeCSC_mvm(p, a, p_out->b);

    FT beta_r = 1.0/beta;
    for (int i = 0; i < p_out->m; i++){
        p_out->b[i] = beta_r * (p->b[i] - p_out->b[i]);
    }

}


