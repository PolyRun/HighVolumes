#include <iostream>
#include "../volume/volume_helper.hpp"


int main(int argc, char **argv){

    CLI cli(argc,argv,"random poly to file");
    CLIFunctionsVolume cliFun(cli);
    
    int n = 12;
    int m = 30;
    int e = 5;
    cliFun.claimOpt('x',"parameters");
    cliFun.add(new CLIF_OptionNumber<int>(&n,'x',"n","12", 2, 100));
    cliFun.add(new CLIF_OptionNumber<int>(&m,'x',"m","30", 4, 1000));
    cliFun.add(new CLIF_OptionNumber<int>(&e,'x',"e","5", 1, 1000));

    
    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    
    std::vector<double> ell(n, 1);
    ell[0] = e;

    
    Polytope *P;
    make_random_poly(ell, m, &P);
    

    /*
    Solved_Body *sb = solved_body_generator()->get("cube_r1.0_3", false);
    Polytope *P = (Polytope *) sb->body[0];
    */    
    //Polytope_print(P);
    // print to vinci format
    std::cout << "random\nH-representation\nbegin\n"
              << P->m << " " << P->n+1 << " real\n";
    for (int i = 0; i < P->m; i++){
        cout << Polytope_get_b(P, i);
        for (int j = 0; j < P->n; j++){
            cout << " " << -Polytope_get_a(P, i, j);
        }
        cout << "\n";
    }
    std::cout << "end\n";

}
