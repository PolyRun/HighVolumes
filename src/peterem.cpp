#include "peterem.hpp"

int main() {
   Polytope* p = Polytope_new(3,6);
   hello(p);

   Polytope* box = Polytope_new_box(4,2);
   hello(box);

   Polytope_free(p);
}

