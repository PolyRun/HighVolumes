#include "volume_helper.hpp"


std::ostream& operator<<(std::ostream& os, const Polytope* p) {
    os << *p; // just refer to print function below.
    return os;
}
std::ostream& operator<<(std::ostream& os, const Polytope& p) {
    os << "Polytope n=" << p.n << ", m=" << p.m << std::endl;

    for(int i=0; i<p.m; i++) {
       // for each constraint

       for(int x=0; x<p.n; x++) {
          os << Polytope_get_a(&p, i, x) << " ";
       }
       os << "| " << Polytope_get_b(&p, i) << std::endl;
    }

    return os;
}


Polytope* Polytope_new_box(int n, int r) {
   Polytope* p = Polytope_new(n, 2*n);

   for(int i=0; i<n; i++) {// for each dim
      Polytope_set_b(p, i,   r);
      Polytope_set_b(p, i+n, r);
      for(int x=0; x<n; x++) {
         Polytope_set_a(p, i,   x, (x==i)?1:0);
         Polytope_set_a(p, i+n, x, (x==i)?-1:0);
      }
   }

   return p;
}


void make_random_poly(const std::vector<double> &ell, int m, Polytope **ret){

  // create n random hyperplanes around ellipsoid ell
  int n = ell.size();
  *ret = Polytope_new(n, m);
  
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution(0,1);

  
  for (int i = 0; i < m; i++){

    // sample uar on unit dim-sphere
    std::vector<double> pnt(n);
    double rad = 0;

    for (int j = 0; j < n; j++) {
      double number = distribution(generator);
      rad += number * number;
      pnt[j] = number;
    }
    rad = sqrt(rad);

    
    // normalize point (on surface) and map point to ellipsoid
    // (note it's no longer uniform at that point but maybe good enough, we could make it uniform by discarding each point with a probability proportional to the distortion factor of the matrix inducing the ellipsoid at that point)
    for (int j = 0; j < n; j++) {
        pnt[j] *= ell[j]/rad;
    }

    // add the tangent plane of ellispoid at pnt to poly
    // note that gradient of the ellipsoid at x is [2*ell_i x_i]_{i = 1 to dim}
    double sum = 0;
    for (int j = 0; j < n; j++){
        Polytope_set_a(*ret, i, j, 2*pnt[j] / (ell[j] * ell[j]));
        sum += Polytope_get_a(*ret, i, j) * pnt[j];
    }
    Polytope_set_b(*ret, i, sum);    
  }

}


