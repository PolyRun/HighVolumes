
// this gives a (mrows x n) to (n x npnts) MMM
int npnts = 4;
int mrows = 4;
int walk_size = 10;
int step_size = 1600;




void PolytopeT_intersect_ref(const void* o, const FT* xs, const FT* ds, FT* t0s, FT* t1s) {
   const PolytopeT* p = (PolytopeT*)o;
   const int n = p->n;
   const int m = p->m;

   for (int i = 0; i < npnts; i++){
       t0s[i] = -FT_MAX;// tmp variables for t0, t1
       t1s[i] = FT_MAX;
   }

   // here we do MMM now
   for(int i = 0; i < m; i++) {
       
       const FT b = PolytopeT_get_b(p, i);

       // TODO: do this in volume_lib_init
       FT *dais = (FT *) calloc(npnts, sizeof(FT));// dotProduct(d,ai,n);
       FT *aixs = (FT *) calloc(npnts, sizeof(FT));// dotProduct(x,ai,n);
       for(int j=0; j<n; j++) {
           FT aij = PolytopeT_get_a(p,i,j);
           for (int p = 0; p < npnts; p++){
               dais[p] += ds[p * n + j] * aij;
               aixs[p] += xs[p * n + j] * aij;
           }
       }
      

       for (int p = 0; p < npnts; p++){
           
           // find intersections y of line with all planes:
           //   y = x + d*t
           //   ai*y = b
           //   
           //   t = (b - ai*x)/(d*ai)
           if (abs(dais[p]) > FT_EPS){
               FT t = (b - aixs[p]) / dais[p];
               if (dais[p] < FT_EPS){
                   t0s[p] = max(t0s[p], t);
               }
               else {
                   t1s[p] = min(t1s[p], t);
               }               
           }
       }
   }
}



void walkCoord_parallel(const int n, const FT rk, int bcount, const void** body, const Body_T** type, FT* x, FT* d, void** cache) {
    const int ws = walk_size; // number of steps for walk
   
    for(int w=0;w<ws;w++) { // take some random steps for x

        // choose npnts random directions
        for (int i = 0; i < npnts; i++){
            int dd = prng_get_random_int_in_range(0,n-1); // pick random dimension
      
        FT t0,t1;
        // ensure do not walk outside of outer ball:
        Ball_intersectCoord(n, rk, x, dd, &t0, &t1);
      
        for(int c=0;c<bcount;c++) {
            FT bt0, bt1;
            type[c]->intersectCoord(body[c], x, dd, &bt0, &bt1, cache[c]);
            t0 = (t0>bt0)?t0:bt0; // max
            t1 = (t1<bt1)?t1:bt1; // min
        }
   
        FT t = prng_get_random_double_in_range(t0,t1);
        x[dd] += t;
        for(int c=0;c<bcount;c++) {
            type[c]->cacheUpdateCoord(body[c], dd, t, cache[c]);
        }

    }
}
