#include <iostream>
#include "../polyvest/vol.h"

extern "C" { // must be included C stlye
#include "../src/poly/volume.h"
#include "../src/poly/cube.h"
#include "../src/ellipsoid.h"
#include "../src/beta_cut.h"
}



void polyvest_convert(Polytope *P, Polyvest_p **Q){

  int n = P->n;
  int m = P->m;
  
  polyvest_p q(m, n);

  for (int i = 0; i < m; i++){
    for (int j = 0; j < n; j++){
      q.matA(Polytope_get_a(P, i, j), i, j);
    }
    q.vecb(Polytope_get_b(P, i), i);
  }

  *Q = q;
}


void test_against_polyvest(Polytope *P){

  int n = P->n;
  
  polyvest_p *Q;
  polyvest_convert(P, &Q);

  
  vec polyvest_ori(n);
  double polyvest_R2;
  Q->genInitE(polyvest_R2, polyvest_ori);

  std::cout << "----------------------------\nPOLYVEST\n";
  std::cout << "initial R2: " << polyvest_R2 << std::endl;
  std::cout << "initial ori:" << std::endl;
  for (int i = 0; i < n; i++){
    cout << polyvest_ori(i) << " ";
  }
  std::cout << endl;
  
  FT R2;
  FT *ori;
  initEllipsoid(P, &R2, &ori);

  
  std::cout << "----------------------------\nHighvolumes\n";
  std::cout << "initial R2: " << R2 << std::endl;
  std::cout << "initial ori:" << std::endl;
  for (int i = 0; i < n; i++){
    cout << ori[i] << " ";
  }
  std::cout << endl;
  
  


}


int main(){
   std::cout << "\n-------------- Test initEllipsoid:\n";

   int n = 5;
   
   Polytope *P, *Q;
   cube(n, &P);

   FT det;

   preprocess(P, &Q, &det);
   

   
   Polytope_free(P);
   Polytope_free(Q);
}
