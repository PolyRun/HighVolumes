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
