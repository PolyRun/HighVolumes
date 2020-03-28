#include <iostream>
#include "../src/poly/volume_helper.hpp"
#include "test_preprocess.hpp"

extern "C" { // must be included C stlye
#include "../src/poly/cube.h"
#include "../src/ellipsoid.h"
#include "../src/beta_cut.h"
}


#define POLYEXP_BASE ((string) "../polyvest/examples/")



void polyvest_convert(Polytope *P, vol::Polyvest_p *Q){

    int n = P->n;
    int m = P->m;  

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            Q->matA(Polytope_get_a(P, i, j), i, j);
        }
        Q->vecb(Polytope_get_b(P, i), i);
    }
  
}


int read_polyvest_p(string filename, Polytope **P){

    ifstream file;
    file.open(filename);


    if (!file.is_open()){
        printf("failed to read polytope");
        return 1;
    }

    int n, m;
    file >> m >> n;

    *P = Polytope_new(n, m);

    FT num;
    for (int i = 0; i < m; i++){
        file >> num;
        Polytope_set_b(*P, i, num);
        for (int j = 0; j < n; j++){
            file >> num;
            Polytope_set_a(*P, i, j, num);
        }
    }

    return 0;
}


void test_init_against_polyvest(Polytope *P){

    int n = P->n;
    int m = P->m;  
  
    vol::Polyvest_p Q(m, n);
    polyvest_convert(P, &Q);

  
    vec polyvest_ori(n);
    double polyvest_R2;
  
    Q.genInitE(polyvest_R2, polyvest_ori);

    std::cout << "---------------- POLYVEST\n";
    std::cout << "initial R2: " << polyvest_R2 << std::endl;
    std::cout << "initial ori:" << std::endl;
    for (int i = 0; i < n; i++){
        cout << polyvest_ori(i) << " ";
    }
    std::cout << endl;
  
    FT R2;
    FT *ori;
    initEllipsoid(P, &R2, &ori);

  
    std::cout << "----------------- Highvolumes\n";
    std::cout << "initial R2: " << R2 << std::endl;
    std::cout << "initial ori:" << std::endl;
    for (int i = 0; i < n; i++){
        cout << ori[i] << " ";
    }
    std::cout << endl;
  
}



void test_preprocess_against_polyvest(Polytope *P){

    int n = P->n;
    int m = P->m;  
  
    vol::Polyvest_p Q(m, n);
    polyvest_convert(P, &Q);

    //Q.A.print();
    //Q.b.print();

    Q.Preprocess();

    Polytope *R;
    FT det;
    preprocess(P, &R, &det);


    std::pair<FT, FT> diff = matrix_diff(R, &Q);

    std::cout << "2-Frobenius of A_P - A_Q:" << std::endl
              << diff.first << std::endl
              << "2-norm of b_P - b_Q:" << std::endl
              << diff.second << std::endl;
             
    
  
}


int main(){
    std::cout << "\n-------------- Test initEllipsoid:\n";

    int n = 5;
   
    //Polytope *P = Polytope_new_box(n, 3);

    Polytope *P;
    string paths[33];
    paths[0] = POLYEXP_BASE + "cc_8_10";
    paths[1] = POLYEXP_BASE + "cc_8_11";
    paths[2] = POLYEXP_BASE + "cross_13";
    paths[3] = POLYEXP_BASE + "cross_7";
    paths[4] = POLYEXP_BASE + "cross_9";
    paths[5] = POLYEXP_BASE + "cube_10";
    paths[6] = POLYEXP_BASE + "cube_10_2";
    paths[7] = POLYEXP_BASE + "cube_14";
    paths[8] = POLYEXP_BASE + "cube_14_2";
    paths[9] = POLYEXP_BASE + "cube_15";
    paths[10] = POLYEXP_BASE + "cube_2";
    paths[11] = POLYEXP_BASE + "cube_20";
    paths[12] = POLYEXP_BASE + "cube_25";
    paths[13] = POLYEXP_BASE + "cube_30";
    paths[14] = POLYEXP_BASE + "cube_35";
    paths[15] = POLYEXP_BASE + "cube_40";
    paths[16] = POLYEXP_BASE + "cube_5";
    paths[17] = POLYEXP_BASE + "cube_80";
    paths[18] = POLYEXP_BASE + "ex_1";
    paths[19] = POLYEXP_BASE + "ex_2";
    paths[20] = POLYEXP_BASE + "fm_6";
    paths[21] = POLYEXP_BASE + "rect_3";
    paths[22] = POLYEXP_BASE + "rh_1";
    paths[23] = POLYEXP_BASE + "rh_2";
    paths[24] = POLYEXP_BASE + "rh_20_40";
    paths[25] = POLYEXP_BASE + "rh_3";
    paths[26] = POLYEXP_BASE + "rh_30_60";
    paths[27] = POLYEXP_BASE + "rh_4";
    paths[28] = POLYEXP_BASE + "rh_40_80";
    paths[29] = POLYEXP_BASE + "simplex_10";
    paths[30] = POLYEXP_BASE + "simplex_14";
    paths[31] = POLYEXP_BASE + "simplex_15";
    paths[32] = POLYEXP_BASE + "simplex_20";


    for (int i = 0; i < 33; i++){

        std::cout << endl << "TESTING " << paths[i] << std::endl;
        
        int err = read_polyvest_p(paths[i], &P);
        if (err){
            return 1;
        }

        FT det;

        //test_init_against_polyvest(P);
        
        test_preprocess_against_polyvest(P);
   
    //preprocess(P, &Q, &det);
   
        Polytope_free(P);
    }
}
