#include <iostream>
#include "../polyvest/vol.h"
#include "../src/poly/volume_helper.hpp"

extern "C" { // must be included C stlye
#include "../src/poly/volume.h"
#include "../src/poly/cube.h"
#include "../src/ellipsoid.h"
#include "../src/beta_cut.h"
}



void polyvest_convert(Polytope *P, vol::Polyvest_p *Q){

    int n = P->n;
    int m = P->m;  

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            Q->matA(Polytope_get_a(P, i, j), i, j);
        }
        Q->vecb(Polytope_get_b(P, i), i);
    }
  
}


void test_init_against_polyvest(Polytope *P){

    int n = P->n;
    int m = P->m;  
  
    vol::Polyvest_p Q(m, n);
    polyvest_convert(P, &Q);

  
    vec polyvest_ori(n);
    double polyvest_R2;
  
    Q.genInitE(polyvest_R2, polyvest_ori);

    std::cout << "---------------- POLYVEST\n";
    std::cout << "initial R2: " << polyvest_R2 << std::endl;
    std::cout << "initial ori:" << std::endl;
    for (int i = 0; i < n; i++){
        cout << polyvest_ori(i) << " ";
    }
    std::cout << endl;
  
    FT R2;
    FT *ori;
    initEllipsoid(P, &R2, &ori);

  
    std::cout << "----------------- Highvolumes\n";
    std::cout << "initial R2: " << R2 << std::endl;
    std::cout << "initial ori:" << std::endl;
    for (int i = 0; i < n; i++){
        cout << ori[i] << " ";
    }
    std::cout << endl;
  
}


void test_preprocess_against_polyvest_output(Polytope *P){

    int n = P->n;
    int m = P->m;  
  
    vol::Polyvest_p Q(m, n);
    polyvest_convert(P, &Q);

    Q.Preprocess();

    Polytope *R;
    FT det;
    preprocess(P, &R, &det);

  
}


int main(){
    std::cout << "\n-------------- Test initEllipsoid:\n";

    int n = 5;
   
    Polytope *P = Polytope_new_box(n, 3);

    FT det;

    //test_init_against_polyvest(P);

    test_preprocess_against_polyvest_output(P);
   
    //preprocess(P, &Q, &det);
   
    Polytope_free(P);
}
