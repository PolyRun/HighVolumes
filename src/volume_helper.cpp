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
