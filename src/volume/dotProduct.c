// imported by volume.c

FT dotProduct_ref(const FT* u, const FT* v, const int n) {
   FT sum = 0.0;
   for(int i=0; i<n; i++) {sum+= u[i]*v[i];}
   return sum;
}

FT dotProduct_vec1(const FT* u, const FT* v, const int n) {
   __m256d sum = _mm256_set1_pd(0.0);
   
   for(int i=0;i<n-3; i+=4) {
      __m256d ui = _mm256_load_pd(u+i);
      __m256d vi = _mm256_load_pd(v+i);
      sum = _mm256_fmadd_pd(ui,vi,sum); // ui*vi + sum
   }
   {// rest
      const int rest = n%4;
      __m256d ui = _mm256_load_pd(u+n-rest);
      __m256d vi = _mm256_load_pd(v+n-rest);
      __m256d s = _mm256_mul_pd(ui,vi);

      __m256d i = _mm256_set_pd(4.0, 3.0, 2.0, 1.0);
      __m256d r = _mm256_set1_pd(rest);
      __m256d c = _mm256_cmp_pd(i,r,_CMP_LE_OS);// i <= r
      //printf("c: %f, %f, %f, %f\n",c[0],c[1],c[2],c[3]);
      s = _mm256_and_pd(c, s);
      //printf("c: %f, %f, %f, %f\n",s[0],s[1],s[2],s[3]);

      sum = _mm256_add_pd(s,sum);
   }

   // sum over sum:
   sum = _mm256_hadd_pd(sum,sum);
   sum = _mm256_permute4x64_pd(sum, 0b00100010);
   sum = _mm256_hadd_pd(sum,sum);
   return sum[0];
}

#include "dotProduct_auto1.c"
