
#ifndef HEADER_BODY_H
#define HEADER_BODY_H

// input body
typedef void (*print_f_t)(const void*);

// input body
typedef void (*free_f_t)(const void*);

// input body
typedef void* (*clone_f_t)(const void*);

// input: body, point x
typedef bool (*inside_f_t)(const void*,const FT*);

// input: body, point x, direction d  -  output t0, t1
// intersections: x+d*t0, x+d*t1
typedef void (*intersect_f_t)(const void*,const FT*,const FT*,FT*,FT*);

// input: body, point x, cooordinate i, cache  -  output t0, t1
typedef void (*intersectCoord_f_t)(const void*,const FT*,const int,FT*,FT*,void*);
typedef FTpair4 (*intersectCoord4_f_t)(const void*,const FT*,const int,void*);
typedef FTpair8 (*intersectCoord8_f_t)(const void*,const FT*,const int,void*);

// input: body - output: number of bites cache required
typedef int (*cacheAlloc_f_t)(const void*);

// input: body, vector x, cache
// for polytopes: computes all m dotproducts Ai x and stores each in cache[i]
typedef void (*cacheReset_f_t)(const void*, const FT*, void*);
typedef void (*cacheReset4_f_t)(const void*, const FT*, void*);
typedef void (*cacheReset8_f_t)(const void*, const FT*, void*);

// input: body, dim d, dx on that dim, cache
typedef void (*cacheUpdateCoord_f_t)(const void*, const int, const FT, void*);
typedef void (*cacheUpdateCoord4_f_t)(const void*, const int, const __m256d, void*);
typedef void (*cacheUpdateCoord8_f_t)(const void*, const int, const FTset8, void*);

// Separation oracle used for preprocessing:
//   If body inside E( (2n)^-2 * A, a):
//       return false
//   else:
//       return true
//       return a plane (v,c) for cutting much of ellipse E(A,a)
//       such that vT * x <= c for all points in body
//       and (2n)^-2 * vT * A * v <= (c - vT * a)^2
//       (cut/touch inner ellipsoid)
//
// input: body, cost/cage ellipsoid
// output: plane (normal v, const c)
//    x in body: vT * x <= c
typedef bool (*shallowCutOracle_f_t)(const void*, const Ellipsoid*, FT*, FT*);

// Transform body after preprocessing
// intput: body_in, body_out, matrix L, vector a, beta.
// old X-space
// new Y-space
// x = (L * y + a)*beta
typedef void (*transform_f_t)(const void*, void*, const Matrix*, const FT*, const FT);

// input: body (ellipsoid or polytope)
// output: radius FT *r and center FT *ori
// compute sphere with center ori and radius r that encloses the body
typedef void (*boundingSphere_f_t)(const void *, FT *, FT *);

// input body, point x
// output: normal n (best effort)
typedef void (*normal_f_t)(const void*, const FT*, FT*);

struct Body_T {
   print_f_t print;
   free_f_t free;
   clone_f_t clone;
   inside_f_t inside;
   intersect_f_t intersect;
   intersectCoord_f_t intersectCoord;
   cacheAlloc_f_t cacheAlloc;
   cacheReset_f_t cacheReset;
   cacheUpdateCoord_f_t cacheUpdateCoord;
   shallowCutOracle_f_t shallowCutOracle;
   transform_f_t transform;
   boundingSphere_f_t boundingSphere;
   normal_f_t normal;
   // 4/8-set:
   intersectCoord4_f_t intersectCoord4;
   cacheReset4_f_t cacheReset4;
   cacheUpdateCoord4_f_t cacheUpdateCoord4;
   intersectCoord8_f_t intersectCoord8;
   cacheReset8_f_t cacheReset8;
   cacheUpdateCoord8_f_t cacheUpdateCoord8;
};


#endif // HEADER_BODY_H
