#include <immintrin.h>
#include <stdbool.h>


bool inside(const double* x) {
   //return x[1] == 1.0;
   if(x[0]*0.1 + x[1]*1.1 + x[6]*4.2 > 3.76) {return false;}
   if(x[4]*0.5 + x[7]*1.2 + x[8]*5.2 > 3.87) {return false;}
   if(x[3]*0.7 + x[9]*1.4 + x[11]*11.2 > 3.47) {return false;}
   return true;
}


void intersect(const double* x, const double* d, double *t0, double *t1) {
   *t0 = *x;
   *t1 = *d;
   return;	
   double t00 = -100000.0;
   double t11 = +100000.0;
   
   // for each:
   
   // 2x dot product: d*a, x*a
   
   const double da = 0.1 * d[0] + 9.7*d[1] + 4.4*d[2];
   const double xa = 0.11 * x[0] + 9.71*x[1] + 4.41*x[2];

   // check if d*a is parallel ==0: jump to next - or make sure t is zero?
   
   const double t = (8.5 - xa) / da;

   t00 += t;
   t11 += t;
   // sub, div
   // cond, min, max

   // end
   *t0 = t00;
   *t1 = t11;
}
