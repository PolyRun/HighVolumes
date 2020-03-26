#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <glpk.h>
#include <chrono>


/*
extern "C" {
#include <stdlib.h>
#include <time.h>
#include "volume.h"
}
*/
using namespace std;

/**
 *\brief this function allow creating a random polytope around an ellipsoid given by ell
 * \param ell an n vector standing for a diagonal matrix
 * \param n the dimension
**/
vector<vector<double>> makepoly(const vector<double> &ell, int n){

  // create n random hyperplanes around ellipsoid ell
  int dim = ell.size();
  vector<vector<double>> poly(n, vector<double> (dim+1));
  
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(seed);
  normal_distribution<double> distribution(0,1);

  
  for (int i = 0; i < n; i++){

    // sample uar on unit dim-sphere
    vector<double> pnt(dim);
    double rad = 0;

    for (int j = 0; j < dim; j++) {
      double number = distribution(generator);
      rad += number * number;
      pnt[j] = number;
    }
    rad = sqrt(rad);

    
    // normalize point (on surface) and map point to ellipsoid
    // (note it's no longer uniform at that point but maybe good enough, we could make it uniform by discarding each point with a probability proportional to the distortion factor of the matrix inducing the ellipsoid at that point)
    for (int j = 0; j < dim; j++) pnt[j] *= ell[j]/rad;

    // add the tangent plane of ellispoid at pnt to poly
    // note that gradient of the ellipsoid at x is [2*ell_i x_i]_{i = 1 to dim}
    double sum = 0;
    for (int j = 0; j < dim; j++){
      poly[i][j] = 2*pnt[j] / (ell[j] * ell[j]);
      sum += poly[i][j] * pnt[j];
    }
    poly[i][dim] = sum;    
  }

  return poly;

}



/**
 * \brief example use of makepoly function reading ellipsoid from ellipse.in
*/
int main(){

  freopen("ellipse.in", "r", stdin);

  // read in ellipsoid (dim x dim matrix) from file
  int dim;
  cin >> dim;
  vector<double> ell(dim);
  for (int i = 0; i < dim; i++){
    cin >> ell[i];
  }


    
  // no idea how many hyperplanes give us good polytopes in expectation:
  // too few and the polytope will likely be ill-conditioned
  // too many and we get basically the ellipsoid back
  int n = 6;

  // TODO: check inequalities are in correct direction!
  vector<vector<double>> poly = makepoly(ell, n); 

  // write polytope, ellipsoid and scaled ellipsoid to out stdout
  for (int i = 0; i < n; i++){
    for (int j = 0; j <= dim; j++){
      cout << poly[i][j] << " ";
    }
    cout << "\n";    
  }

  /* write plotting data for testing */
      
}
