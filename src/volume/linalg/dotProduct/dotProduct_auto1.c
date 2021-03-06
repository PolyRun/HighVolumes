// Code Generated by dotProduct_gen.py
#include "dotProduct.h"
FT dotProduct_auto1(const FT* u, const FT* v, const int n) {
   switch(n) {
      case 0: {
         return 0;
         break;
      }//case 0
      case 1: {
         // m = 1
         // rest = 1
         const __m256i mask = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
         const __m256d u0 = _mm256_maskload_pd(u+0,mask);
         const __m256d v0 = _mm256_maskload_pd(v+0,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 1
      case 2: {
         // m = 1
         // rest = 2
         const __m256i mask = _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u0 = _mm256_maskload_pd(u+0,mask);
         const __m256d v0 = _mm256_maskload_pd(v+0,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 2
      case 3: {
         // m = 1
         // rest = 3
         const __m256i mask = _mm256_set_epi64x(0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u0 = _mm256_maskload_pd(u+0,mask);
         const __m256d v0 = _mm256_maskload_pd(v+0,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 3
      case 4: {
         // m = 1
         // rest = 0
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         __m256d sum = _mm256_mul_pd(u0,v0);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 4
      case 5: {
         // m = 2
         // rest = 1
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
         const __m256d u1 = _mm256_maskload_pd(u+4,mask);
         const __m256d v1 = _mm256_maskload_pd(v+4,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         sum = _mm256_fmadd_pd(u1,v1,sum);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 5
      case 6: {
         // m = 2
         // rest = 2
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u1 = _mm256_maskload_pd(u+4,mask);
         const __m256d v1 = _mm256_maskload_pd(v+4,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         sum = _mm256_fmadd_pd(u1,v1,sum);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 6
      case 7: {
         // m = 2
         // rest = 3
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256i mask = _mm256_set_epi64x(0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u1 = _mm256_maskload_pd(u+4,mask);
         const __m256d v1 = _mm256_maskload_pd(v+4,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         sum = _mm256_fmadd_pd(u1,v1,sum);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 7
      case 8: {
         // m = 2
         // rest = 0
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         __m256d sum = _mm256_mul_pd(u0,v0);
         sum = _mm256_fmadd_pd(u1,v1,sum);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 8
      case 9: {
         // m = 3
         // rest = 1
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
         const __m256d u2 = _mm256_maskload_pd(u+8,mask);
         const __m256d v2 = _mm256_maskload_pd(v+8,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 9
      case 10: {
         // m = 3
         // rest = 2
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u2 = _mm256_maskload_pd(u+8,mask);
         const __m256d v2 = _mm256_maskload_pd(v+8,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 10
      case 11: {
         // m = 3
         // rest = 3
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256i mask = _mm256_set_epi64x(0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u2 = _mm256_maskload_pd(u+8,mask);
         const __m256d v2 = _mm256_maskload_pd(v+8,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 11
      case 12: {
         // m = 3
         // rest = 0
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 12
      case 13: {
         // m = 4
         // rest = 1
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
         const __m256d u3 = _mm256_maskload_pd(u+12,mask);
         const __m256d v3 = _mm256_maskload_pd(v+12,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 13
      case 14: {
         // m = 4
         // rest = 2
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u3 = _mm256_maskload_pd(u+12,mask);
         const __m256d v3 = _mm256_maskload_pd(v+12,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 14
      case 15: {
         // m = 4
         // rest = 3
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256i mask = _mm256_set_epi64x(0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u3 = _mm256_maskload_pd(u+12,mask);
         const __m256d v3 = _mm256_maskload_pd(v+12,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 15
      case 16: {
         // m = 4
         // rest = 0
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 16
      case 17: {
         // m = 5
         // rest = 1
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
         const __m256d u4 = _mm256_maskload_pd(u+16,mask);
         const __m256d v4 = _mm256_maskload_pd(v+16,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 17
      case 18: {
         // m = 5
         // rest = 2
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u4 = _mm256_maskload_pd(u+16,mask);
         const __m256d v4 = _mm256_maskload_pd(v+16,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 18
      case 19: {
         // m = 5
         // rest = 3
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256i mask = _mm256_set_epi64x(0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u4 = _mm256_maskload_pd(u+16,mask);
         const __m256d v4 = _mm256_maskload_pd(v+16,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 19
      case 20: {
         // m = 5
         // rest = 0
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 20
      case 21: {
         // m = 6
         // rest = 1
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
         const __m256d u5 = _mm256_maskload_pd(u+20,mask);
         const __m256d v5 = _mm256_maskload_pd(v+20,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 21
      case 22: {
         // m = 6
         // rest = 2
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u5 = _mm256_maskload_pd(u+20,mask);
         const __m256d v5 = _mm256_maskload_pd(v+20,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 22
      case 23: {
         // m = 6
         // rest = 3
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256i mask = _mm256_set_epi64x(0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u5 = _mm256_maskload_pd(u+20,mask);
         const __m256d v5 = _mm256_maskload_pd(v+20,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 23
      case 24: {
         // m = 6
         // rest = 0
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256d u5 = _mm256_load_pd(u+20);
         const __m256d v5 = _mm256_load_pd(v+20);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 24
      case 25: {
         // m = 7
         // rest = 1
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256d u5 = _mm256_load_pd(u+20);
         const __m256d v5 = _mm256_load_pd(v+20);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
         const __m256d u6 = _mm256_maskload_pd(u+24,mask);
         const __m256d v6 = _mm256_maskload_pd(v+24,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_fmadd_pd(u6,v6,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 25
      case 26: {
         // m = 7
         // rest = 2
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256d u5 = _mm256_load_pd(u+20);
         const __m256d v5 = _mm256_load_pd(v+20);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u6 = _mm256_maskload_pd(u+24,mask);
         const __m256d v6 = _mm256_maskload_pd(v+24,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_fmadd_pd(u6,v6,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 26
      case 27: {
         // m = 7
         // rest = 3
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256d u5 = _mm256_load_pd(u+20);
         const __m256d v5 = _mm256_load_pd(v+20);
         const __m256i mask = _mm256_set_epi64x(0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u6 = _mm256_maskload_pd(u+24,mask);
         const __m256d v6 = _mm256_maskload_pd(v+24,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_fmadd_pd(u6,v6,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 27
      case 28: {
         // m = 7
         // rest = 0
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256d u5 = _mm256_load_pd(u+20);
         const __m256d v5 = _mm256_load_pd(v+20);
         const __m256d u6 = _mm256_load_pd(u+24);
         const __m256d v6 = _mm256_load_pd(v+24);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_fmadd_pd(u6,v6,sum);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 28
      case 29: {
         // m = 8
         // rest = 1
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256d u5 = _mm256_load_pd(u+20);
         const __m256d v5 = _mm256_load_pd(v+20);
         const __m256d u6 = _mm256_load_pd(u+24);
         const __m256d v6 = _mm256_load_pd(v+24);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFF);
         const __m256d u7 = _mm256_maskload_pd(u+28,mask);
         const __m256d v7 = _mm256_maskload_pd(v+28,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_fmadd_pd(u6,v6,sum);
         sum2 = _mm256_fmadd_pd(u7,v7,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 29
      case 30: {
         // m = 8
         // rest = 2
         const __m256d u0 = _mm256_load_pd(u+0);
         const __m256d v0 = _mm256_load_pd(v+0);
         const __m256d u1 = _mm256_load_pd(u+4);
         const __m256d v1 = _mm256_load_pd(v+4);
         const __m256d u2 = _mm256_load_pd(u+8);
         const __m256d v2 = _mm256_load_pd(v+8);
         const __m256d u3 = _mm256_load_pd(u+12);
         const __m256d v3 = _mm256_load_pd(v+12);
         const __m256d u4 = _mm256_load_pd(u+16);
         const __m256d v4 = _mm256_load_pd(v+16);
         const __m256d u5 = _mm256_load_pd(u+20);
         const __m256d v5 = _mm256_load_pd(v+20);
         const __m256d u6 = _mm256_load_pd(u+24);
         const __m256d v6 = _mm256_load_pd(v+24);
         const __m256i mask = _mm256_set_epi64x(0, 0, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
         const __m256d u7 = _mm256_maskload_pd(u+28,mask);
         const __m256d v7 = _mm256_maskload_pd(v+28,mask);
         __m256d sum = _mm256_mul_pd(u0,v0);
         __m256d sum2 = _mm256_mul_pd(u1,v1);
         sum = _mm256_fmadd_pd(u2,v2,sum);
         sum2 = _mm256_fmadd_pd(u3,v3,sum2);
         sum = _mm256_fmadd_pd(u4,v4,sum);
         sum2 = _mm256_fmadd_pd(u5,v5,sum2);
         sum = _mm256_fmadd_pd(u6,v6,sum);
         sum2 = _mm256_fmadd_pd(u7,v7,sum2);
         sum = _mm256_add_pd(sum,sum2);
         sum = _mm256_hadd_pd(sum,sum);
         return sum[0]+sum[2];
         break;
      }//case 30
      default:{
         FT sum = 0.0;
         for(int i=0; i<n; i++) {sum+= u[i]*v[i];}
         return sum;
         break;
      }//case default
   }//switch
}//dotProduct_auto1
