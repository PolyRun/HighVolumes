#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <glpk.h>
#include <chrono>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Gmpq.h>

typedef CGAL::Gmpq IT;
typedef CGAL::Gmpq ET;
typedef CGAL::Quadratic_program<IT> Prog;
typedef CGAL::Quadratic_program_solution<ET> Sol;


/*
extern "C" {
#include <stdlib.h>
#include <time.h>
#include "volume.h"
}
*/
using namespace std;

/*
bool checkfeasibility(glp_prob *lp, const vector<double> &ell, double scale){
  
  int n = ell.size();
  
  for (int i = 0; i < n; i++){
    double scaled_coef = ell[i] * scale;
    double inf_s_c = 1./scaled_coef;
    glp_set_obj_coef(lp, i+1, inf_s_c);
  }


  //disable msg output
  glp_smcp parm;
  glp_init_smcp(&parm);
  parm.msg_lev = GLP_MSG_ERR;
  
  glp_simplex(lp, &parm);
    
  double val = glp_get_obj_val(lp);
  
  if (val <= 1) { // all points in polytope are inside ellipsoid
    return true;
  }
  else {
    return false;
  }
}
*/

bool checkfeasibility(Prog *qp, const vector<double> &ell, double scale){
  
  int n = ell.size();
  
  for (int i = 0; i < n; i++){
    double scaled_coef = ell[i] * scale;
    double inf_s_c = 1./(scaled_coef * scaled_coef);
    qp->set_d(i, i, -inf_s_c);
  }

  Sol s = CGAL::solve_quadratic_program(*qp, ET());
  assert (s.solves_quadratic_program(*qp));

  cout << s;
  
  if (s.is_infeasible() || s.is_unbounded())
    throw("polytope degenerate");
    
  auto val = s.objective_value();
  
  if (val >= -1) { // all points in polytope are inside ellipsoid
    return true;
  }
  else {
    return false;
  }

  
}


double scaleell(const vector<vector<double>> &poly, const vector<double> &ell){

  int m = poly.size(), n = poly[0].size()-1;
  bool contains = false;
  double lo = 1, hi = -1;
  double scale = lo;

  // TODO: how to set eps appropriately? ideal scaling will be in [res - eps, res], this might have impact on runtime in high dimensions
  double eps = .0001;
  
  Prog qp(CGAL::SMALLER, false, 0, false, 0);

  for (int i = 0; i < m; i++){
    for (int j = 0; j < n; j++){
      qp.set_a(j, i, poly[i][j]);
    }
    qp.set_b(i, poly[i][n]);
  }
  
  // exponential search scaling factor such that scale*ell contains poly
  while (!contains){
    scale *= 2;
    contains = checkfeasibility(&qp, ell, scale);
  }
  hi = scale;

  // binary search for best value between lo and hi
  while (hi - lo > eps){
    scale = lo + (hi - lo)/2;
    contains = checkfeasibility(&qp, ell, scale);

    if (contains) { // all points in polytope are inside ellipsoid
      hi = scale;
    }
    else { 
      lo = scale;
    }
  }
  
  return hi;



}



/*
double scaleell(const vector<vector<double>> &poly, const vector<double> &ell){


  int m = poly.size(), n = poly[0].size()-1;
  bool contains = false;
  double lo = 1, hi = -1;
  double scale = lo;


  // TODO: how to set eps appropriately? ideal scaling will be in [res - eps, res], this might have impact on runtime in high dimensions
  double eps = .0001;

  //init GLPK
  glp_prob *lp;
  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX); // maximize objective
  glp_add_rows(lp, m);
  glp_add_cols(lp, n);


  // load constraints
  int *ind = new int[n + 1];
  double *val = new double[n + 1];
  for (int i = 0; i < m; i++){
    for (int j = 0; j < n; j++){
      ind[j+1] = j+1;
      val[j+1] = poly[i][j];
    }
    glp_set_mat_row(lp, i+1, n, ind, val);
    glp_set_row_bnds(lp, i+1, GLP_UP, 0, poly[i][n]);
  }
  delete []ind, delete []val;
  
  for (int i = 0; i < n; i++) {
    glp_set_col_bnds(lp, i+1, GLP_FR, 0, 0); // all free variables
  }



  checkfeasibility(lp, ell, 1);
  
  // exponential search scaling factor such that scale*ell contains poly
  while (!contains){
    scale *= 2;
    contains = checkfeasibility(lp, ell, scale);
  }
  hi = scale;

  // binary search for best value between lo and hi
  while (hi - lo > eps){
    scale = lo + (hi - lo)/2;
    contains = checkfeasibility(lp, ell, scale);

    if (contains) { // all points in polytope are inside ellipsoid
      hi = scale;
    }
    else { 
      lo = scale;
    }
  }
  
  glp_delete_prob(lp);
  return hi;

}
*/


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



/*
  create random polytope around an ellipsoid
  @returns: writes polytope and circumscribing ellispoid to file
  @parameter: the in file should contain an n doubles d_i, then A := diag(d_i) induces the ellispoid (transform unit ball)
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

  double scale = scaleell(poly, ell);

  cout << scale << "\n";
  
  freopen("out.out", "w", stdout);

  /*

  // write polytope, ellipsoid and scaled ellipsoid to out stdout

  for (int i = 0; i < n; i++){
    for (int j = 0; j <= dim; j++){
      cout << poly[i][j] << " ";
    }
    cout << "\n";    
  }

  cout << "\n";
  for (int i = 0; i < dim; i++){
    cout << ell[i] << " ";
  }
  cout << "\n\n";
  for (int i = 0; i < dim; i++){
    cout << scale*ell[i] << " ";
  }
  cout << "\n";
  */

  /* write plotting data for testing */


  vector<double> ell_scaled(dim);
  for (int i = 0; i < dim; i++) ell_scaled[i] = ell[i]*scale;
  
  int npnts = 1000;
  double xmax = 1./(ell[0] * ell[0]), xmin = -xmax,
    xmax_scaled = 1./(ell_scaled[0] * ell_scaled[0]), xmin_scaled = -xmax_scaled,
    x_start = -3, x_end = 3;

  for (int i = 0; i < npnts; i++){
    double x = x_start + i*(x_end - x_start)/npnts;
    
    cout << x;

    // output the polytope
    double minabove = 1 << 30, maxbelow = -1 << 30;
    for (int i = 0; i < n; i++){
      double y = (poly[i][dim] - poly[i][0] * x)/ poly[i][1];
      if (poly[i][1] > 0 && y < minabove)
	minabove = y;
      else if (poly[i][1] < 0 && y > maxbelow)
	maxbelow = y;
    }

    if (minabove >= maxbelow){
      cout << ", " << minabove << ", " << maxbelow;
    }
    else {
      cout << ", nan, nan";
    }

    // output the ellipse at x
    //    if (x >= xmin && x <= xmax){
      double y = ell[1]*ell[1] * (1 - (x*x)/(ell[0]*ell[0]));
      cout << ", " << sqrt(y) << ", " << -sqrt(y);
      //}

    // output the scaled ellipse at x
    //if (x >= xmin_scaled && x <= xmax_scaled){
      double y_scaled = ell_scaled[1]*ell_scaled[1] * (1 - (x*x)/(ell_scaled[0]*ell_scaled[0]));
      cout << ", " << sqrt(y_scaled) << ", " << -sqrt(y_scaled);
      //}

    cout << "\n";
  }
      
}
