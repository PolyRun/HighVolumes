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
    
    for (int i = p->col_start[d]; i < p->col_start[d+1] && p->row_idx[i] > -1; i++){
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



void PolytopeCSC_intersectCoord_cached_withb(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache){

    // Note: in this case the cache stores cache[i] = b[i] - Ai * x
    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    FT *b_Aix = (FT *) cache;
    
    FT t00 = -FT_MAX;
    FT t11 = FT_MAX;
    
    for (int i = p->col_start[d]; i < p->col_start[d+1] && p->row_idx[i] > -1; i++){
        FT dai = p->A[i];

        FT t = (b_Aix[p->row_idx[i]]) / dai;

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




void PolytopeCSC_intersectCoord_cached_vec(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache){

    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    FT *b_Aix = (FT *) cache;

    int colstrt = p->col_start[d];
    FT *A = &p->A[colstrt];
    int *row = &p->row_idx[colstrt];
    
    // elems is a multiple of 4!
    // thus we need to consider exactly c := |col|/4 vectors
    int elems = p->col_start[d+1] - colstrt;

    __m128i rowmask = _mm_set1_epi32(0xffffffff);
    
    __m256d zeros = _mm256_set1_pd(0.0);

    __m256d t00 = _mm256_set1_pd(-FT_MAX);
    __m256d t11 = _mm256_set1_pd(FT_MAX);

    // the first c-1 vectors have all valid elements
    int i = 0;
    for (; i < elems - 7; i+=4) {

        __m128i rows = _mm_maskload_epi32(row + i, rowmask);
        
        __m256d b_Aix_vec = _mm256_i32gather_pd(b_Aix, rows, 8);
        __m256d dai = _mm256_load_pd(A + i);

        __m256d t = _mm256_div_pd(b_Aix_vec, dai);

        // <comp>_OS -> signal nans
        // n_mask[i] is 0xff..ff if t[i] < 0 and 0 else
        __m256d n_mask = _mm256_cmp_pd(dai, zeros, _CMP_LT_OS);

        __m256d t00_tmp = _mm256_blendv_pd(t00, t, n_mask);
        __m256d t11_tmp = _mm256_blendv_pd(t, t11, n_mask);

        t00 = _mm256_max_pd(t00, t00_tmp);
        t11 = _mm256_min_pd(t11, t11_tmp);
    }

    // only last vector has invalid entries
    __m128i rows = _mm_maskload_epi32(row + i, rowmask);
    //__m128i rows = _mm_load_si128((__m128i const*) (row + i));
    // note that rows is negative iff the entry is not valid, thus we can use it as mask (if we negate it first)
    __m256d mask = _mm256_cvtepi32_pd(rows);
    // negate mask, is there a better way?
    //mask = _mm256_fmsub_pd(zeros, zeros, mask);
        
    __m256d b_Aix_vec = _mm256_i32gather_pd(b_Aix, rows, 8);
    __m256d dai = _mm256_load_pd(A + i);

    __m256d t = _mm256_div_pd(b_Aix_vec, dai);

    // <comp>_OS -> signal nans
    // n_mask[i] is 0xff..ff if t[i] < 0 and 0 else
    // p_mask[i] is 0xff..ff if t[i] >= 0 and 0 else
    // is there a smarter way to invert masks?
    // note: logic is cheaper than cmp and blend!
    __m256d n_mask = _mm256_cmp_pd(t, zeros, _CMP_LT_OS);
    __m256d p_mask = _mm256_cmp_pd(t, zeros, _CMP_GE_OS);
    __m256d nv_mask = _mm256_andnot_pd(mask, n_mask);
    __m256d pv_mask = _mm256_andnot_pd(mask, p_mask);
    
    __m256d t00_tmp = _mm256_blendv_pd(t00, t, nv_mask);
    //__m256d t00_tmptmp = _mm256_blendv_pd(t00_tmp, t00, mask);
    __m256d t11_tmp = _mm256_blendv_pd(t11, t, pv_mask);
    //__m256d t11_tmptmp = _mm256_blendv_pd(t11_tmp, t11, mask);

    t00 = _mm256_max_pd(t00, t00_tmp);
    t11 = _mm256_min_pd(t11, t11_tmp);

    FT t0_tmp = t00[0];
    FT t0_tmptmp = (t00[1] > t0_tmp) ? t00[1] : t0_tmp;
    FT t0_tmptmptmp = (t00[2] > t0_tmptmp) ? t00[2] : t0_tmptmp;
    *t0 = (t00[3] > t0_tmptmptmp) ? t00[3] : t0_tmptmptmp;

    FT t1_tmp = t11[0];
    FT t1_tmptmp = (t11[1] < t1_tmp) ? t11[1] : t1_tmp;
    FT t1_tmptmptmp = (t11[2] < t1_tmptmp) ? t11[2] : t1_tmptmp;
    *t1 = (t11[3] < t1_tmptmptmp) ? t11[3] : t1_tmptmptmp;

}


static inline void loop(const int *rowi, const FT *b_Aix, const FT *Ai, const __m128i rowmask, const __m256d zeros, __m256d *t00, __m256d *t11){
    
    __m128i rows = _mm_maskload_epi32(rowi, rowmask);
        
    __m256d b_Aix_vec = _mm256_i32gather_pd(b_Aix, rows, 8);
    __m256d dai = _mm256_load_pd(Ai);

    __m256d t = _mm256_div_pd(b_Aix_vec, dai);

    // <comp>_OS -> signal nans
    // n_mask[i] is 0xff..ff if t[i] < 0 and 0 else
    __m256d n_mask = _mm256_cmp_pd(dai, zeros, _CMP_LT_OS);

    __m256d t00_tmp = _mm256_blendv_pd(*t00, t, n_mask);
    __m256d t11_tmp = _mm256_blendv_pd(t, *t11, n_mask);

    *t00 = _mm256_max_pd(*t00, t00_tmp);
    *t11 = _mm256_min_pd(*t11, t11_tmp);
}


void PolytopeCSC_intersectCoord_cached_vec_inline(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache){

    const PolytopeCSC *p = (PolytopeCSC *) o;
    int n = p->n;
    int m = p->m;
    FT *b_Aix = (FT *) cache;

    int colstrt = p->col_start[d];
    FT *A = &p->A[colstrt];
    int *row = &p->row_idx[colstrt];
    
    // elems is a multiple of 4!
    // thus we need to consider exactly c := |col|/4 vectors
    int elems = p->col_start[d+1] - colstrt;

    __m128i rowmask = _mm_set1_epi32(0xffffffff);
    
    __m256d zeros = _mm256_set1_pd(0.0);

    __m256d t00 = _mm256_set1_pd(-FT_MAX);
    __m256d t11 = _mm256_set1_pd(FT_MAX);

    // the first c-1 vectors have all valid elements
    int i = 0;
    for (; i < elems - 11; i+=4) {
        
        loop(row+i, b_Aix, A+i, rowmask, zeros, &t00, &t11);
        loop(row+i+4, b_Aix, A+i+4, rowmask, zeros, &t00, &t11);
    }

    if (i < elems - 7){
        loop(row+i, b_Aix, A+i, rowmask, zeros, &t00, &t11);
        i+=4;
    }

    // only last vector has invalid entries
    __m128i rows = _mm_maskload_epi32(row + i, rowmask);
    //__m128i rows = _mm_load_si128((__m128i const*) (row + i));
    // note that rows is negative iff the entry is not valid, thus we can use it as mask (if we negate it first)
    __m256d mask = _mm256_cvtepi32_pd(rows);
    // negate mask, is there a better way?
    //mask = _mm256_fmsub_pd(zeros, zeros, mask);
        
    __m256d b_Aix_vec = _mm256_i32gather_pd(b_Aix, rows, 8);
    __m256d dai = _mm256_load_pd(A + i);

    __m256d t = _mm256_div_pd(b_Aix_vec, dai);

    // <comp>_OS -> signal nans
    // n_mask[i] is 0xff..ff if t[i] < 0 and 0 else
    // p_mask[i] is 0xff..ff if t[i] >= 0 and 0 else
    // is there a smarter way to invert masks?
    // note: logic is cheaper than cmp and blend!
    __m256d n_mask = _mm256_cmp_pd(t, zeros, _CMP_LT_OS);
    __m256d p_mask = _mm256_cmp_pd(t, zeros, _CMP_GE_OS);
    __m256d nv_mask = _mm256_andnot_pd(mask, n_mask);
    __m256d pv_mask = _mm256_andnot_pd(mask, p_mask);
    
    __m256d t00_tmp = _mm256_blendv_pd(t00, t, nv_mask);
    //__m256d t00_tmptmp = _mm256_blendv_pd(t00_tmp, t00, mask);
    __m256d t11_tmp = _mm256_blendv_pd(t11, t, pv_mask);
    //__m256d t11_tmptmp = _mm256_blendv_pd(t11_tmp, t11, mask);

    t00 = _mm256_max_pd(t00, t00_tmp);
    t11 = _mm256_min_pd(t11, t11_tmp);

    FT t0_tmp = t00[0];
    FT t0_tmptmp = (t00[1] > t0_tmp) ? t00[1] : t0_tmp;
    FT t0_tmptmptmp = (t00[2] > t0_tmptmp) ? t00[2] : t0_tmptmp;
    *t0 = (t00[3] > t0_tmptmptmp) ? t00[3] : t0_tmptmptmp;

    FT t1_tmp = t11[0];
    FT t1_tmptmp = (t11[1] < t1_tmp) ? t11[1] : t1_tmp;
    FT t1_tmptmptmp = (t11[2] < t1_tmptmp) ? t11[2] : t1_tmptmp;
    *t1 = (t11[3] < t1_tmptmptmp) ? t11[3] : t1_tmptmptmp;

}
