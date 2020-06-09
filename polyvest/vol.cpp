/***********************************************************************
 *  This code is part of PolyVest.
 *
 *  Copyright (C) 2013, 2016 Cunjing Ge, Institute of Software, Chinese 
 *  Academy of Sciences, Beijing, China. All rights reserved. 
 *  E-mail: <gecj@ios.ac.cn>.
 *
 *  PolyVest is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  PolyVest is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PolyVest. If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#include "vol.h"
#include "glpk.h"
#include <iomanip>

#define PI 3.1415926536
#define FLOATWIDTH 15
//#define DEBUG_MSG
//#define DEBUG
//#define PRINT_T
//#define PRINT_TMI
using namespace vol;

//double abs(double x){
//	return (x > 0) ? x : -x;
//}

/**
   compute volume of B(0, n)
**/
double uballVol(int n){
    double vol = 1;
    if (n % 2 == 1){
        int k = (n - 1) / 2;
        vol *= pow(2, n);
        for (int i = 1; i < k + 1; i++) vol *= PI * i;
        for (int i = 1; i < n + 1; i++) vol /= i;
    }else{
        int k = n / 2;
        for (int i = 1; i < k + 1; i++) vol *= PI / i;
    }
    return vol;
}

/*********** Delete Redundent Hyperplanes ***********/
void Polyvest_p::checkHPs(){
    if (check_planes_off) return;

    //init GLPK
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, m);
    glp_add_cols(lp, n);

    //disable msg output
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR;

    //load constraints
    int *ind = new int[n + 1];
    double *val = new double[n + 1];
    for (int i = 1; i < m + 1; i++){
        for (int j = 1; j < n + 1; j++){
            ind[j] = j, val[j] = A(i - 1, j - 1);
        }
        glp_set_mat_row(lp, i, n, ind, val);
        glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
    }
    delete []ind, delete []val;
    for (int i = 1; i < n + 1; i++)
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);

    //feasiblity check
    int num[2];
    for (int i = 1; i < m + 1;){
        glp_set_row_bnds(lp, i, GLP_LO, b(i - 1) + 0.00001, 0);
        glp_set_obj_coef(lp, 1, 1);
        for (int j = 1; j < n; j++)
            glp_set_obj_coef(lp, j + 1, 0);
        glp_simplex(lp, &parm);
        if (glp_get_status(lp) == GLP_NOFEAS){
            cout << "\n\nHYPERPLANE REMOVED IN POLYVEST\n"
                 << "expect different performance\n\n";
            num[1] = i;
            glp_del_rows(lp, 1, num);
            A.shed_row(i - 1);
            b.shed_row(i - 1);
            m--;
        }else{
            glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
            i++;
        }
    }
    cout << "Hyperplanes Left: " << m << endl;
    glp_delete_prob(lp);
}


void Polyvest_p::genInitE(double &R2, vec &Ori){
    R2 = 0, Ori.zeros();

    //init GLPK
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, m);
    glp_add_cols(lp, n);


    //disable msg output
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR;

    //load constraints
    int *ind = new int[n + 1];
    double *val = new double[n + 1];
    
    for (int i = 1; i < m + 1; i++){
        for (int j = 1; j < n+1; j++){
            ind[j] = j, val[j] = A(i - 1, j - 1);
        }
        glp_set_mat_row(lp, i, n, ind, val);
        glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
    }
    delete []ind, delete []val;
    for (int i = 1; i < n + 1; i++)
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);

    //get bounds
    for (int i = 0; i < n; i++){
        double max, min;
        for (int j = 0; j < n; j++)
            glp_set_obj_coef(lp, j + 1, 0);

        glp_set_obj_coef(lp, i + 1, 1);
        glp_simplex(lp, &parm);
        max = glp_get_obj_val(lp);
        for (int j = 0; j < n; j++)
            Ori(j) += glp_get_col_prim(lp, j + 1);

        glp_set_obj_coef(lp, i + 1, -1);
        glp_simplex(lp, &parm);
        min = -glp_get_obj_val(lp);
        for (int j = 0; j < n; j++)
            Ori(j) += glp_get_col_prim(lp, j + 1);

        R2 += (max - min) * (max - min);
    }
    Ori = Ori / (2 * n);
	
    glp_delete_prob(lp);
}



/**
   - remove redundant hyperplanes using checkHPs
   - compute minimum enclosing ellipsoid using genInitE
   - transform body given the axes-lengths of ellipsoid
   - store determinant of transformation for volume scaling
**/
void Polyvest_p::Preprocess_hacked(){

    std::cout << std::fixed << setprecision(FLOATWIDTH);
    
    //checkHPs();

    double c1 = (2 * pow(n, 2) + pow(1 - n / beta_r, 2)) * (1 - 1.0 / pow(beta_r, 2)) / (2 * pow(n, 2) - 2);
    //double c1 = pow(n, 2) * (1 - 1.0 / pow(beta_r, 2)) / (pow(n, 2) - 1);
    double c2 = (1 - n / beta_r) / (n + 1);
    double c3 = beta_r * beta_r;
    double c4 = 2 * c2 / (1 - 1.0 / beta_r);


    
#ifdef DEBUG
    cout << "c1: " << c1 << "\n" << endl;
    cout << "c2: " << c2 << "\n" << endl;
    cout << "c3: " << c3 << "\n" << endl;
    cout << "c4: " << c4 << "\n" << endl;
#endif
    
    //init E(R2I, 0), T = R2I, ori = 0.
    mat T;
    vec ori(n);
    double R2;
    genInitE(R2, ori);
    T.eye(n, n);
    T = R2 * T;


#ifdef DEBUG_MSG
        std::cout << "--------------- POLYVEST\n"
                  << "First ellipsoid approximation\n"
                  << "T:\n";
        T.raw_print();
        std::cout << "\ncenter:\n";
        ori.raw_print();
#endif

    
    vec distance = zeros<vec>(m);
    vec tm = zeros<vec>(m);

    int counter = 0;
    while (++counter > 0){
        int i;
		
        //check if ori in polytope
        distance = b - A * ori;

#ifdef DEBUG
        //cout << "ROUND " << counter << endl;
        //distance.t().raw_print();
#endif

        for (i = 0; i < m; i++)
            if (distance(i) < 0){
#ifdef DEBUG
                //cout << i << " in LOOP 1, distance[" << i << "] = " << distance(i) << " < 0" << endl;
#endif
                tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
                //cout << "tm" << i << ": " << tm(i) << endl;
                break;
            }

        
        if (i == m){
            //check if small ellipsoid contained in polytope
            for (i = 0; i < m; i++){
#ifdef DEBUG
                //cout << i << " in LOOP 2" << endl;
#endif
                tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
                if (c3 * distance(i) * distance(i) - tm(i) < 0) {
#ifdef DEBUG
                    //cout << c3 << "*" << distance(i) << "*" << distance(i) << "-" << tm(i) << "=" << c3*distance(i)*distance(i)-tm(i) << "< 0" << endl;
#endif
                    break;
                }
                else {
#ifdef DEBUG
                    //cout << c3 << "*" << distance(i) << "*" << distance(i) << "-" << tm(i) << "=" << c3*distance(i)*distance(i)-tm(i) << ">= 0" << endl;
#endif
                }
            }
        }
        

        
#ifdef PRINT_TMI
        //printf("ROUND %d\n", counter);
        tm.t().raw_print();
#endif


        //terminate if E satisfies two criterions
        if (i == m) break;
		
        vec t = T * A.row(i).t() / sqrt(tm(i));
        ori = ori - t * c2;
        T = c1 * (T - c4 * t * t.t());

#ifdef PRINT_T
        ori.raw_print();
        cout << endl;
        T.raw_print();
        cout << endl;
#endif

    }

    
	

	
    //apply affine transformation
    //mat Trans = chol(T);
	
    mat Trans;
    try{
        Trans = chol(T);
    }catch (const std::runtime_error& ex){
        cout << "The input polytope is degenerated or non-existed and the volume is 0." << endl;
        exit(1);		
    }

    
    //b = beta_r * (b - A * ori);
    //A = A * Trans.t();


    rowvec exp(n);
    exp.ones();
    for (int i = 0; i < n; i++){
        B[i] = b / A.col(i);
        Ai[i] = A / (A.col(i) * exp);
    }
	
    determinant = 1;//det(Trans) / pow(beta_r, n);


#ifdef DEBUG_MSG
    
        std::cout << "Final ellipsoid\n"
                  << "T:\n";
        T.raw_print();
        
        std::cout << "\ncenter:\n";
        ori.raw_print();
        
        std::cout << "\nTrans:\n";
        Trans.raw_print();
        
        std::cout << "Transformed Poly:\n"
                  << "A:\n";
        A.raw_print();

        std::cout << "b:\n";
        b.raw_print();
        
        std::cout << "Determinant:\n"
                  << determinant << std::endl;
        std::cout << "The number of iterations: " << counter << std::endl
                  << "^^^^^^^^^^^^^^^^^ END POLYVEST" << std::endl;
#endif

}




double Polyvest_p::EstimateVol(int coef = 1600){
    int k, i, maxk, maxl = 0;

    const long stepsz = coef * l; //size of sampling

    double *alpha = new double[l];
    long *volK = new long[l];
    memset(alpha, 0, l * sizeof(double));
    memset(volK, 0, l * sizeof(long));
	
    x.zeros();
    int steps = 0;
    for (k = l - 2; k >= 0; k--){
        for (i = volK[k + 1]; i < stepsz; i++){
            steps++;
            double m = walk(k);
            if (m < r2[0]) volK[0]++;
            else if (m < r2[k])
                volK[(int)trunc(n * log(m) / (log((double)2) * 2)) + 1]++;
        }
        for (i = 0; i < k; i++){
            volK[k] += volK[i];
        }
        if (volK[k] < stepsz){
            alpha[k] = (double)(stepsz) / volK[k];
            x = x / pow((double)2, (double)1 / n);
        }else alpha[k] = 1;
    }
    cout << "Polyvest did " << steps << " steps\n";
    cout << "There are " << l << " shells\n";
    vol = uballVol(n) * determinant;

    for (i = 0; alpha[i] > 1 && i < l - 1; i++){
        vol *= alpha[i];
    }

    delete []alpha;
    delete []volK;

    return vol;
}

// Original
void Polyvest_p::Preprocess(){

    checkHPs();

    double c1 = (2 * pow(n, 2) + pow(1 - n / beta_r, 2)) * (1 - 1.0 / pow(beta_r, 2)) / (2 * pow(n, 2) - 2);
    //double c1 = pow(n, 2) * (1 - 1.0 / pow(beta_r, 2)) / (pow(n, 2) - 1);
    double c2 = (1 - n / beta_r) / (n + 1);
    double c3 = beta_r * beta_r;
    double c4 = 2 * c2 / (1 - 1.0 / beta_r);

    cout << "find init ellipsoid...\n";
    //init E(R2I, 0), T = R2I, ori = 0.
    mat T;
    vec ori(n);
    double R2;
    genInitE(R2, ori);
    T.eye(n, n);
    T = R2 * T;

    vec distance = zeros<vec>(m);
    vec tm = zeros<vec>(m);
    
    cout << "start stepping...\n";

    int counter = 0;
    while (++counter > 0){
        int i;
		
        //check if ori in polytope
        distance = b - A * ori;
        for (i = 0; i < m; i++) {
            if (distance(i) < 0){
                tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
                break;
            }
        }
        
        if (i == m){
            //check if small ellipsoid contained in polytope
            for (i = 0; i < m; i++){
                tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
                if (c3 * distance(i) * distance(i) - tm(i) < 0){
                    break;
                }
            }
        }
		
        //terminate if E satisfies two criterions
        if (i == m) break;
		
        vec t = T * A.row(i).t() / sqrt(tm(i));
        ori = ori - t * c2;
        T = c1 * (T - c4 * t * t.t());
    }
    cout << "Steps: " << counter << "\n";
	
    if (!msg_off){ 
        cout << "R^2: " << R2 << endl << "Origin: " << endl;
        ori.print();
    }
	
    //apply affine transformation
    //mat Trans = chol(T);
	
    mat Trans;
    try{
        Trans = chol(T);
    }catch (const std::runtime_error& ex){
        cout << "The input polytope is degenerated or non-existed and the volume is 0." << endl;
        exit(1);		
    }
	
    //cout << Trans << endl;
    b = beta_r * (b - A * ori);
    A = A * Trans.t();

    if (!msg_off) cout << "The number of iterations: " << counter << endl;

    rowvec exp(n);
    exp.ones();
    for (int i = 0; i < n; i++){
        B[i] = b / A.col(i);
        Ai[i] = A / (A.col(i) * exp);
    }
	
    determinant = det(Trans) / pow(beta_r, n);
}



/**
   - generate one random point on layer k, i.e. for comparing K_k with K_{k-1}
   - coordinate walk! they choose direction dir as randi!!   
**/
double Polyvest_p::walk(int k){
    double r, max, min, C = 0;
    int dir = randi(n);
    vec::col_iterator it, itA, end;

    // find distance to boundary of polytope in direction dir and -dir
    for (it = x.begin(), end = x.end(); it != end; ++it) C += *it * *it;
    C -= x(dir) * x(dir);
    r = sqrt(r2[k + 1] - C);
    max = r - x(dir), min = -r - x(dir);

    //A(x + t v) <= b
    //Av t <= b - Ax

	
    vec bound = B[dir] - Ai[dir] * x;
    for (it = bound.begin(), end = bound.end(), itA = A.begin_col(dir); it != end; ++it, ++itA){
        if (*itA > 0){
            if (*it < max) max = *it;
        }else if (*itA < 0)
            if (*it > min) min = *it; 
    }

    // choose random point on line given by dir
    double t = x(dir) + randd(max - min) + min;
    x(dir) = t;

    return (C + t * t);
}
