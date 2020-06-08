#include "linalg.h"

Ball_intersect_f_t Ball_intersect = Ball_intersect_ref;
void Ball_intersect_ref(const int n, const FT r, const FT* x, const FT* d, FT* t0, FT* t1) {
   // y = x + d*t
   // y^2 = r^2
   //
   // x^2 - r^2 + 2x*d*t + d^2*t^2 = 0
   
   FT x2 = squaredNorm(x,n);
   FT d2 = squaredNorm(d,n);
   FT xd = dotProduct(x,d,n);

   FT a = d2;
   FT ainv = 1.0 / a;
   FT b = 2.0*xd;
   FT c = x2 - r*r;

   FT detSqrt = sqrt(b*b - 4.0*a*c);
   
   *t1 = (-b + detSqrt) * 0.5 * ainv;
   *t0 = (-b - detSqrt) * 0.5 * ainv;
}


Ball_intersectCoord_f_t Ball_intersectCoord = Ball_intersectCoord_ref;
void Ball_intersectCoord_ref(const int n, const FT r, const FT* x, const int d, FT* t0, FT* t1) {
   
   FT x2 = squaredNorm(x,n);
   const FT d2 = 1.0;
   FT xd = x[d]; // dot product with unit vector dim d

   const FT a = d2;
   const FT ainv = 1.0 / a;
   FT b = 2.0*xd;
   FT c = x2 - r*r;

   FT detSqrt = sqrt(b*b - 4.0*a*c);
   
   *t1 = (-b + detSqrt) * 0.5 * ainv;
   *t0 = (-b - detSqrt) * 0.5 * ainv;
}

FT Ball_volume(const int n, const FT r) {
   FT rn = pow(r,n);
   
   int nhalf = n / 2;
   bool odd = (n % 2) == 1;
   
   FT pinhalf = pow(M_PI,nhalf);

   FT i = nhalf;
   if(odd) {i+=0.5;}
   
   FT fact = 1.0;
   while(i>0) {fact*=i; i-=1.0;}
   
   return pinhalf * rn / fact;
}

int ceil_cache(const int n, const int b) {
   return (n*b + 31) / 32 * 32;
}


squaredNorm_cached_f_t squaredNorm_cached = squaredNorm_cached_ref;
squaredNorm_cached_reset_f_t squaredNorm_cached_reset = squaredNorm_cached_reset_ref;
squaredNorm_cached_update_f_t squaredNorm_cached_update = squaredNorm_cached_update_ref;


FT squaredNorm_cached_ref(const FT* v, const int n, const FT* cache) {
   // not cached
   return squaredNorm(v,n);
}
void squaredNorm_cached_reset_ref(const FT* v, const int n, FT* cache) {
   // not cached
}
void squaredNorm_cached_update_ref(const FT* v, const int d, const FT dx, const int n, FT* cache) {
   // not cached
}

FT squaredNorm_cached_refc(const FT* v, const int n, const FT* cache) {
   return *cache;
}
void squaredNorm_cached_reset_refc(const FT* v, const int n, FT* cache) {
   *cache = squaredNorm(v,n);
}
void squaredNorm_cached_update_refc(const FT* v, const int d, const FT dx, const int n, FT* cache) {
   FT c = *cache;
   FT x = v[d];
   FT y = x - dx;
   *cache = c - y*y + x*x; // subtract old square, add new square
}

Ball_intersectCoord_cached_f_t Ball_intersectCoord_cached = Ball_intersectCoord_cached_ref;

FTpair Ball_intersectCoord_cached_ref(const int n, const FT r, const FT* x,const int d, FT* cache) {
   FT x2 = squaredNorm_cached(x,n,cache);
   const FT d2 = 1.0;
   FT xd = x[d]; // dot product with unit vector dim d

   const FT a = d2;
   const FT ainv = 1.0 / a;
   FT b = 2.0*xd;
   FT c = x2 - r*r;

   FT detSqrt = sqrt(b*b - 4.0*a*c);
   
   FTpair tp;
   tp.t1 = (-b + detSqrt) * 0.5 * ainv;
   tp.t0 = (-b - detSqrt) * 0.5 * ainv;
   return tp;
}


// ----------------------------------4/8 sets:

__m256d squaredNorm_cached4(const FT* v, const int n, const FT* cache) {
   return _mm256_load_pd(cache);
}
void squaredNorm_cached4_reset(const FT* v, const int n, FT* cache) {
   __m256d acc = _mm256_set1_pd(0.0);
   // TODO: make more efficient!
   for(int i=0;i<n;i++) {
      __m256d vi = _mm256_load_pd(v+4*i);
      acc = _mm256_fmadd_pd(vi,vi,acc);
   }
   _mm256_store_pd(cache,acc);
}
void squaredNorm_cached4_update(const FT* v, const int d, const __m256d dx, const int n, FT* cache) {
   __m256d c = _mm256_load_pd(cache);
   __m256d x = _mm256_load_pd(v+4*d);
   __m256d y = _mm256_sub_pd(x,dx);
   __m256d y2 = _mm256_mul_pd(y,y);
   c = _mm256_fmadd_pd(x,x,c);
   c = _mm256_sub_pd(c,y2);
   _mm256_store_pd(cache,c);
}

FTset8 squaredNorm_cached8(const FT* v, const int n, const FT* cache) {
   __m256d c0 = _mm256_load_pd(cache);
   __m256d c1 = _mm256_load_pd(cache+4);
   FTset8 res = {c0,c1};
   return res;
}
void squaredNorm_cached8_reset(const FT* v, const int n, FT* cache) {
   __m256d acc0 = _mm256_set1_pd(0.0);
   __m256d acc1 = _mm256_set1_pd(0.0);
   // TODO: make more efficient!
   for(int i=0;i<n;i++) {
      __m256d vi0 = _mm256_load_pd(v+8*i);
      __m256d vi1 = _mm256_load_pd(v+8*i+4);
      acc0 = _mm256_fmadd_pd(vi0,vi0,acc0);
      acc1 = _mm256_fmadd_pd(vi1,vi1,acc1);
   }
   _mm256_store_pd(cache,acc0);
   _mm256_store_pd(cache+4,acc1);
}  

void squaredNorm_cached8_update(const FT* v, const int d, const FTset8 dx, const int n, FT* cache) {
   __m256d c0 = _mm256_load_pd(cache);
   __m256d c1 = _mm256_load_pd(cache+4);

   __m256d x0 = _mm256_load_pd(v+8*d);
   __m256d x1 = _mm256_load_pd(v+8*d+4);

   __m256d y0 = _mm256_sub_pd(x0,dx.set0);
   __m256d y1 = _mm256_sub_pd(x1,dx.set1);
   
   __m256d y20 = _mm256_mul_pd(y0,y0);
   __m256d y21 = _mm256_mul_pd(y1,y1);
   
   c0 = _mm256_fmadd_pd(x0,x0,c0);
   c1 = _mm256_fmadd_pd(x1,x1,c1);
   
   c0 = _mm256_sub_pd(c0,y20);
   c1 = _mm256_sub_pd(c1,y21);
   
   _mm256_store_pd(cache,c0);
   _mm256_store_pd(cache+4,c1);
}


FTpair4 Ball_intersectCoord_cached4(const int n, const FT r, const FT* x,const int d, FT* cache) {
   //FT x2 = squaredNorm_cached(x,n,cache);
   __m256d x2 = squaredNorm_cached4(x,n,cache);
   
   //FT xd = x[d]; // dot product with unit vector dim d
   __m256d xd = _mm256_load_pd(x+4*d); // dot product with unit vector dim d
   
   __m256d two = _mm256_set1_pd(2.0);
   __m256d four = _mm256_set1_pd(4.0);
   __m256d mhalf = _mm256_set1_pd(-0.5);
   __m256d rr = _mm256_set1_pd(r);

   //FT b = 2.0*xd;
   __m256d b = _mm256_mul_pd(two,xd);
   //FT c = x2 - r*r;
   __m256d c = _mm256_fmsub_pd(rr,rr,x2);// c negative!
   
   //FT detSqrt = sqrt(b*b - 4.0*a*c);
   __m256d tmp = _mm256_mul_pd(four,c);
   __m256d tmp_ = _mm256_fmadd_pd(b,b,tmp);
   __m256d detSqrt = _mm256_sqrt_pd(tmp_);
   
   //tp.t1 = (-b + detSqrt) * 0.5 * ainv;
   //tp.t0 = (-b - detSqrt) * 0.5 * ainv;
   __m256d hi = _mm256_sub_pd(b,detSqrt);
   __m256d lo = _mm256_add_pd(b,detSqrt);
   FTpair4 tp;
   tp.hi0 = _mm256_mul_pd(hi,mhalf);
   tp.low0  = _mm256_mul_pd(lo,mhalf);
   return tp;
}
FTpair8 Ball_intersectCoord_cached8(const int n, const FT r, const FT* x,const int d, FT* cache) {
   FTset8 xx22 = squaredNorm_cached8(x,n,cache);
   __m256d x20 = xx22.set0;
   __m256d x21 = xx22.set1;

   __m256d xd0 = _mm256_load_pd(x+8*d); // dot product with unit vector dim d
   __m256d xd1 = _mm256_load_pd(x+8*d+4); // dot product with unit vector dim d
   
   __m256d two = _mm256_set1_pd(2.0);
   __m256d four = _mm256_set1_pd(4.0);
   __m256d mhalf = _mm256_set1_pd(-0.5);
   __m256d rr = _mm256_set1_pd(r);

   __m256d b0 = _mm256_mul_pd(two,xd0);
   __m256d b1 = _mm256_mul_pd(two,xd1);
   __m256d c0 = _mm256_fmsub_pd(rr,rr,x20);// c negative!
   __m256d c1 = _mm256_fmsub_pd(rr,rr,x21);// c negative!
   
   __m256d tmp0 = _mm256_mul_pd(four,c0);
   __m256d tmp1 = _mm256_mul_pd(four,c1);
   __m256d tmp_0 = _mm256_fmadd_pd(b0,b0,tmp0);
   __m256d tmp_1 = _mm256_fmadd_pd(b1,b1,tmp1);
   __m256d detSqrt0 = _mm256_sqrt_pd(tmp_0);
   __m256d detSqrt1 = _mm256_sqrt_pd(tmp_1);
   
   __m256d hi0 = _mm256_sub_pd(b0,detSqrt0);
   __m256d hi1 = _mm256_sub_pd(b1,detSqrt1);
   __m256d lo0 = _mm256_add_pd(b0,detSqrt0);
   __m256d lo1 = _mm256_add_pd(b1,detSqrt1);
   FTpair8 tp;
   tp.hi0 = _mm256_mul_pd(hi0,mhalf);
   tp.hi1 = _mm256_mul_pd(hi1,mhalf);
   tp.low0  = _mm256_mul_pd(lo0,mhalf);
   tp.low1  = _mm256_mul_pd(lo1,mhalf);
   return tp;
}

ArbitraryExpNum ArbitraryExpNum_new(FT val) {
   ArbitraryExpNum a;
   a.num = val;
   a.numLog = log(val);
   return a;
}
ArbitraryExpNum ArbitraryExpNum_mul(ArbitraryExpNum a, FT val) {
   ArbitraryExpNum b;
   b.num = a.num * val;
   b.numLog = a.numLog + log(val);
   return b;
}
ArbitraryExpNum ArbitraryExpNum_mul2(ArbitraryExpNum a, ArbitraryExpNum b) {
   ArbitraryExpNum c;
   c.num = a.num * b.num;
   c.numLog = a.numLog + b.numLog;
   return c;
}


void ArbitraryExpNum_print(ArbitraryExpNum a) {
   FT b10 = a.numLog / log(10);
   int exp = floor(b10);
   FT m = pow(10,b10-exp);
   if(exp>300 || exp<-300 || a.num > (double)1e300 || a.num==0) {
      printf("[m: %lf, e: %d]",m, exp);
   }else{
      printf("[%.10e, log: %.10e, m: %lf, e: %d]",a.num,a.numLog, m, exp);
   }
}


// default, doesn't use shell_cache -> does nothing
void shell_cache_init_nocache(FT *cache, FT r0, int l, FT stepFac){}
shell_cache_init_f_t shell_cache_init = shell_cache_init_nocache; 

void shell_cache_init_ref(FT *cache, FT r0, int l, FT stepFac){
    cache[0] = r0;
    // store squaredNorms of shell radii
    for (int i = 1; i <= l; i++){
        cache[i] = cache[i-1]/stepFac;
        cache[i-1] *= cache[i-1];
    }
    cache[l] *= cache[l];
}


int shell_idx_nocache(FT x2, FT r0, FT stepFac, int l, FT *cache){
    const FT mmm = log(x2/(r0*r0))*0.5/(-log(stepFac));
    const int mm = ceil(mmm);
    return (mm>0)?mm:0; // find index of balls
}
shell_idx_f_t shell_idx = shell_idx_nocache;

int shell_idx_ref(FT x2, FT r0, FT stepFac, int l, FT *cache){
    // do linear search in shell_cache for first index larger than x2
    int i = 0;
    while(cache[i] < x2){ i++; }
    return i;
       
}

int shell_idx_binary(FT x2, FT r0, FT stepFac, int l, FT *cache){
    // do binary search in shell_cache for first index larger than x2
    int hi = l;
    int lo = 0;

    while (hi - lo > 0) {
        int mid = lo + (hi - lo)/2;
        if (cache[mid] < x2){
            lo = mid+1;
        }
        else {
            hi = mid;
        }
    }
    return lo;
}

