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
        .boundingSphere = PolytopeJIT_bounding_ref
};

PolytopeJIT* PolytopeJIT_new(int n, int m) {
   PolytopeJIT* o = (PolytopeJIT*) malloc(sizeof(PolytopeJIT));
   o->n = n;
   o->m = m;
   o->inside = NULL;
   o->intersect = NULL;
   o->intersectCoord = NULL;
   return o;
}

PolytopeJIT *Polytope_to_PolytopeJIT(const Polytope *p) {
   const int n = p->n;
   const int m = p->m;
   PolytopeJIT* o = PolytopeJIT_new(n,m);
   PolytopeJIT_print(o);
   
   PolytopeJIT_generate_inside_ref(p,o);
   PolytopeJIT_generate_intersect_ref(p,o);
   PolytopeJIT_generate_intersectCoord_ref(p,o);
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
   printf("intersectCoord: %p\n",p->intersectCoord);
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
   return p->intersect(x,d,t0,t1);
}

void PolytopeJIT_intersectCoord_ref(const void* o, const FT* x, const int d, FT* t0, FT* t1, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   const int n = p->n;
   const int m = p->m;
   
   assert(false && "intersectCoord not implemented for PolytopeJIT");
}

int  PolytopeJIT_cacheAlloc_ref(const void* o) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   // allocate one FT per inequality: store dot product: dot(ai,x)
   return p->m * sizeof(FT);
}
void PolytopeJIT_cacheReset_ref(const void* o, const FT* x, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(false && "cacheReset not implemented for PolytopeJIT");
}

void PolytopeJIT_cacheUpdateCoord_ref(const void* o, const int d, const FT dx, void* cache) {
   const PolytopeJIT* p = (PolytopeJIT*)o;
   assert(false && "cacheUpdateCoord not implemented for PolytopeJIT");
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

