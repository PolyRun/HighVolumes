#include <immintrin.h>
#include <stdbool.h>


bool inside(const double* x) {
   //return x[1] == 1.0;
   if(x[0]*0.1 + x[1]*1.1 + x[6]*4.2 > 3.76) {return false;}
   if(x[4]*0.5 + x[7]*1.2 + x[8]*5.2 > 3.87) {return false;}
   if(x[3]*0.7 + x[9]*1.4 + x[11]*11.2 > 3.47) {return false;}
   return true;
}



