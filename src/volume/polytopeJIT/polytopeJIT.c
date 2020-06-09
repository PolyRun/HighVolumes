#include "polytopeJIT.h"

Body_T PolytopeJIT_T = {
	.print = PolytopeJIT_print,
        .free = PolytopeJIT_free,
        .clone = PolytopeJIT_clone,
        .inside = PolytopeJIT_inside_ref,
        .intersect = PolytopeJIT_intersect_ref,
        .intersectCoord = PolytopeJIT_intersectCoord_ref,
	.cacheAlloc = PolytopeJIT_cacheAlloc_ref,
	.cacheReset = PolytopeJIT_cacheReset_ref,
	.cacheUpdateCoord = PolytopeJIT_cacheUpdateCoord_ref,
        .shallowCutOracle = PolytopeJIT_shallowCutOracle_ref,
	.transform = PolytopeJIT_transform_ref,
        .boundingSphere = PolytopeJIT_bounding_ref,
        .intersectCoord4 = PolytopeJIT_intersectCoord4_ref,
        .intersectCoord8 = PolytopeJIT_intersectCoord8_ref,
	.cacheReset4 = PolytopeJIT_cacheReset4_ref,
	.cacheReset8 = PolytopeJIT_cacheReset8_ref,
	.cacheUpdateCoord4 = PolytopeJIT_cacheUpdateCoord4_ref,
	.cacheUpdateCoord8 = PolytopeJIT_cacheUpdateCoord8_ref,
};

PolytopeJIT_Generator PolytopeJIT_generator = pjit_single_rax;

PolytopeJIT* PolytopeJIT_new(int n, int m) {
   PolytopeJIT* o = (PolytopeJIT*) malloc(sizeof(PolytopeJIT));
   o->n = n;
   o->m = m;
   o->inside = NULL;
   o->intersect = NULL;
   o->cacheReset = NULL;
   o->intersectCoord = NULL;
   o->cacheUpdateCoord = NULL;
   o->intersectCoord_bytes = 0;
   o->cacheUpdateCoord_bytes = 0;
   
   o->cacheReset4 = NULL;
   o->cacheReset8 = NULL;
   o->intersectCoord4 = NULL;
   o->intersectCoord8 = NULL;
   o->cacheUpdateCoord4 = NULL;
   o->cacheUpdateCoord8 = NULL;
   
   return o;
}

PolytopeJIT *Polytope_to_PolytopeJIT(const Polytope *p) {
   const int n = p->n;
   const int m = p->m;
   
   PolytopeJIT* o = PolytopeJIT_new(n,m);

   // count nzA
   o->nzA = 0;
   for(int i=0;i<p->n;i++) {
      for(int j=0;j<p->m;j++) {
         FT aij = Polytope_get_a(p,j,i);
	 if(aij != 0.0) {
	    o->nzA++;
	 }
      }
   }

   //PolytopeJIT_print(o);
   
   PolytopeJIT_generate_inside_ref(p,o);
   PolytopeJIT_generate_intersect_ref(p,o);
   PolytopeJIT_generate_cacheReset_ref(p,o);
   PolytopeJIT_generate_intersectCoord_ref(p,o);
   PolytopeJIT_generate_cacheUpdateCoord_ref(p,o);
   return o;
}

void PolytopeJIT_free(const void* o) {
   PolytopeJIT* p = (PolytopeJIT*)o;
   free(p);
}

void* PolytopeJIT_clone(const void* o) {
   PolytopeJIT* old = (PolytopeJIT*)o;
   assert(false && "Cannot clone PolytopeJIT!");
   return NULL;
}

void PolytopeJIT_print(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   printf("PolytopeJIT: n=%d, m=%d\n",p->n,p->m);

   printf("inside: %p\n",p->inside);
   printf("intersect: %p\n",p->intersect);
   printf("cacheReset: %p\n",p->cacheReset);
   printf("intersectCoord: %p   (%ld bytes)\n",p->intersectCoord, p->intersectCoord_bytes);
   printf("cacheUpdateCoord: %p   (%ld bytes)\n",p->cacheUpdateCoord, p->cacheUpdateCoord_bytes);
   printf("nonZero entries in A: %d (%f)\n",p->nzA,(double)p->nzA/(p->n*p->m));
}

bool PolytopeJIT_inside_ref(const void* o, const FT* v) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(p->inside && "inside function must be generated PolytopeJIT");
   return p->inside(v);
}


void PolytopeJIT_intersect_ref(const void* o, const FT* x, const FT* d, FT* t0, FT* t1) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   const int m = p->m;
   
   assert(p->intersect && "intersect function must be generated PolytopeJIT");
   p->intersect(x,d,t0,t1);
   return;
}

void PolytopeJIT_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   const int m = p->m;
   assert(p->intersectCoord && "intersectCoord function must be generated PolytopeJIT");
   //printf("intersectCoord %d\n",d);
   //for(int i=0;i<n;i++) {printf(" %f",x[i]);} printf(" x\n");
   //for(int i=0;i<m;i++) {printf(" %f",((double*)cache)[i]);} printf(" c\n");
   p->intersectCoord(d,t0,t1,cache);
   //printf("-> %lf %lf\n",*t0,*t1);
}

FTpair4 PolytopeJIT_intersectCoord4_ref(const void* o, const FT* x, const int d, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   const int m = p->m;
   assert(p->intersectCoord4 && "intersectCoord4 function must be generated PolytopeJIT");
   
   FTpair4 tp = {{0,0,0,0},{0,0,0,0}};
   //printf("intersectCoord4\n");
   p->intersectCoord4(d,&tp.low0,&tp.hi0,cache);
   //printf("%lf %lf %lf %lf - low0\n",tp.low0[0],tp.low0[1],tp.low0[2],tp.low0[3]);
   //printf("%lf %lf %lf %lf - hi0\n",tp.hi0[0],tp.hi0[1],tp.hi0[2],tp.hi0[3]);
   return tp;
}
FTpair8 PolytopeJIT_intersectCoord8_ref(const void* o, const FT* x, const int d, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   const int m = p->m;
   assert(p->intersectCoord8 && "intersectCoord8 function must be generated PolytopeJIT");
   
   FTpair8 tp = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
   //printf("intersectCoord8\n");
   p->intersectCoord8(d,&tp.low0,&tp.hi0,cache);
   //printf("%lf %lf %lf %lf - low0\n",tp.low0[0],tp.low0[1],tp.low0[2],tp.low0[3]);
   //printf("%lf %lf %lf %lf - low1\n",tp.low1[0],tp.low1[1],tp.low1[2],tp.low1[3]);
   //printf("%lf %lf %lf %lf - hi0\n",tp.hi0[0],tp.hi0[1],tp.hi0[2],tp.hi0[3]);
   //printf("%lf %lf %lf %lf - hi1\n",tp.hi1[0],tp.hi1[1],tp.hi1[2],tp.hi1[3]);
   return tp;
}
void PolytopeJIT_cacheReset4_ref(const void* o, const FT* x, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(p->cacheReset4 && "cacheReset4 function must be generated PolytopeJIT");

   //for(int i=0;i<4*p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" before\n");
   
   p->cacheReset4(x,cache);
   
   //for(int i=0;i<4*p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" after\n");
}
void PolytopeJIT_cacheReset8_ref(const void* o, const FT* x, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(p->cacheReset8 && "cacheReset8 function must be generated PolytopeJIT");
   //for(int i=0;i<8*p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" before\n");
   p->cacheReset8(x,cache);
   //for(int i=0;i<8*p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" after\n");
}

void PolytopeJIT_cacheUpdateCoord4_ref(const void* o, const int d, const __m256d dx, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(p->cacheUpdateCoord4 && "cacheReset4 function must be generated PolytopeJIT");
   
   //printf("cacheUpdateCoord4 %d %f\n",d,dx);
   //for(int i=0;i<4*p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" before\n");
   p->cacheUpdateCoord4(d,dx,cache);
   //for(int i=0;i<4*p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" after\n");
}

void PolytopeJIT_cacheUpdateCoord8_ref(const void* o, const int d, const FTset8 dx, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(p->cacheUpdateCoord8 && "cacheReset8 function must be generated PolytopeJIT");
   
   //printf("cacheUpdateCoord %d %f\n",d,dx);
   //for(int i=0;i<p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" before\n");
   p->cacheUpdateCoord8(d,dx.set0,dx.set1,cache);
   //for(int i=0;i<p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" after\n");
}

int  PolytopeJIT_cacheAlloc_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   // allocate one FT per inequality: store dot product: dot(ai,x)
   return (p->m+4) * sizeof(FT);
   // +4 just in case something writes over end...
}
void PolytopeJIT_cacheReset_ref(const void* o, const FT* x, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(p->cacheReset && "cacheReset function must be generated PolytopeJIT");
   p->cacheReset(x,cache);
}

void PolytopeJIT_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(p->cacheUpdateCoord && "cacheUpdateCoord function must be generated PolytopeJIT");
   
   //printf("cacheUpdateCoord %d %f\n",d,dx);
   //for(int i=0;i<p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" before\n");
   p->cacheUpdateCoord(d,dx,cache);
   //for(int i=0;i<p->m;i++) {printf(" %f",((double*)cache)[i]);}
   //printf(" after\n");
}

bool PolytopeJIT_shallowCutOracle_ref(const void* o, const Ellipsoid* e, FT* v, FT* c) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   const int m = p->m;
   
   assert(false && "shallowCutOracle not implemented for PolytopeJIT");
   
   // no half-plane violated inner ellipse
   return false;
}

void PolytopeJIT_transform_ref(const void* o_in, void* o_out, const Matrix* L, const FT* a, const FT beta) {
   const PolytopeJIT* p_in = (PolytopeJIT*)o_in;
   PolytopeJIT* p_out = (PolytopeJIT*)o_out;
   const int n = p_in->n;
   const int m = p_in->m;
   
   assert(false && "transform not implemented for PolytopeJIT");
}

void PolytopeJIT_bounding_ref(const void *B, FT *R2, FT *ori) {
   assert(false && "bounding not implemented for PolytopeJIT");
}

