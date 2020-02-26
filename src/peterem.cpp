#include "peterem.hpp"

int main() {
   Polytope* p = Polytope_new(3,6);

   hello(p);

   Polytope_free(p);
}

