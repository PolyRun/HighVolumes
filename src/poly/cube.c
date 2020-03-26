#include "cube.h"

void cube(int dim, Polytope **P){
    int n = dim;
    int m = 2*dim;
    
    *P = Polytope_new(n,m);
    for (int i = 0; i < n; i++){
        // x_i >= 0
        Polytope_set_a(*P, 2*i, i, -1);
        Polytope_set_b(*P, 2*i, 0);
        //(*P)->data[(2*i)*(n+1) + i] = -1;
        //(*P)->data[(2*i)*(n+1) + n] = 0;
        // x_i <= 1
        Polytope_set_a(*P, 2*i+1, i, 1);
        Polytope_set_b(*P, 2*i+1, 1);
        //(*P)->data[(2*i+1)*(n+1) + i] = 1;
        //(*P)->data[(2*i+1)*(n+1) + n] = 1;
    }
}
