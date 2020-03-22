#include "volume.h"

Polytope* Polytope_new(int n, int m) {
   Polytope* o = (Polytope*) malloc(sizeof(Polytope));
   o->data = (FT*) malloc(sizeof(FT)*(n+1)*m);
   o->n = n;
   o->m = m;

   return o;
}


void Polytope_free(Polytope* p) {
   free(p->data);
   free(p);
}


void Polytope_print(Polytope *p){
    int n = p->n;
    int m = p->m;
    printf("Printing poly n: %d, m: %d at %p\n\n", n, m, p);
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            printf("%f ", p->data[i*(n+1)+j]);
        }
        printf("<= %f\n", p->data[i*(n+1)+n]);
    }
    printf("\n");
}
