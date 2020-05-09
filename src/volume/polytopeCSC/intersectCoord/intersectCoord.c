#include "intersectCoord.h"
#include <immintrin.h>


void PolytopeCSC_intersectCoord_cached_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache){

    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    FT *Aix = (FT *) cache;
    FT *b = p->b;
    
    FT t00 = -FT_MAX;
    FT t11 = FT_MAX;
    
    for (int i = p->col_start[d]; i < p->col_start[d+1], p->row_idx[i] > -1; i++){
        FT bi = b[p->row_idx[i]];
        FT dai = p->A[i];

        FT t = (bi - Aix[p->row_idx[i]]) / dai;

        if (dai < 0.0){
            t00 = (t00 > t) ? t00 : t;
        }
        else {
            t11 = (t11 < t) ? t11 : t;
        }
    }

    *t0 = t00;
    *t1 = t11;
}



/*
void PolytopeCSC_intersectCoord_cached_vec(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache){

    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    FT *Aix = (FT *) cache;
    FT *b = p->b;
    
    FT t00 = -FT_MAX;
    FT t11 = FT_MAX;

    int colstrt = col_start[d];
    FT *A = &p->A[colstrt];
    // elems is a multiple of 4!
    int elems = col_start[d+1] - colstrt;

    _m256d srcb = _mm256_set1_pd(0.0);
    _m256d srcdai = _mm256_set1_pd(-1.0);
    
    for (int i = 0; i < elems - 3; i+=4) {

        __m128i rows = _mm128_load_epi32(p->row_idx[d] + i);
        // note that rows is negative iff the entry is not valid, thus we can use it as mask (negative doubles have first bit set)
        // no, need to invert it first!
        __m256d mask = _mm256_cvtepi32_pd(rows); 
        
        __m256d bi = _mm256_mask_i32gather_pd(srcb, b, rows, mask, 8);
        __m256d Aix = _mm256_mask_i32gather_pd(
        __m256d dai = _mm256_load_pd(A + i);

        
        FT bi = b[p->row_idx[i]];
        FT dai = p->A[i];

        FT t = (bi - Aix[p->row_idx[i]]) / dai;

        if (dai < 0.0){
            t00 = (t00 > t) ? t00 : t;
        }
        else {
            t11 = (t11 < t) ? t11 : t;
        }

        i+=4;
    }

    *t0 = t00;
    *t1 = t11;
}
*/
