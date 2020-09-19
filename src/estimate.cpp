#include "volume/volume_helper.hpp"



int main(int argc, char *argv[]){
    CLI_LONG cli(argc,argv,"main");
    CLI_LONG_Functions cliFun(cli);
    initVolumeFunctions(cliFun);
    
    int r = 100;
    int warmup = 0;
    double time_ci_alpha;
    double results_ci_alpha;
    cliFun.add_long(new CLIF_OptionNumber<int>(&r,'r',"r","100", 1, 100000));
    cliFun.add_long(new CLIF_OptionNumber<int>(&warmup,'w',"warmup","0", 0, 100000));
    cliFun.add_long(new CLIF_OptionNumber<double>(&time_ci_alpha,0,"time_ci_alpha","0.95", 0, 1));
    cliFun.add_long(new CLIF_OptionNumber<double>(&results_ci_alpha,0,"results_ci_alpha","0.95", 0, 1));

    std::string filename = "";
    cliFun.add_long(new CLIF_MandatoryString(&filename, 'f', "filename"));

    
    bool polytopeOptimize = false;
    cliFun.add_long(new CLIF_Option<bool>(&polytopeOptimize,'o',"polytopeOptimize","false", {
                                                     {"false",{false,"-"}},
						     {"true", {true, "Sort constraints to optimize access pattern"}} }));

    bool printBody = false;
    cliFun.add_long(new CLIF_Option<bool>(&printBody,'b',"printBody","false", {
                                                     {"false",{false,"-"}},
						     {"true", {true, "Print body before benchmark is run."}} }));

    int polytopeType = 0;
    cliFun.add_long(new CLIF_Option<int>(&polytopeType,'t',"polytopeType","0",
                                    {
                                     {"0",{0, "Polytope format / rows"}},
                                     {"1",{1, "PolytopeT format / columns"}},
                                     {"2",{2, "PolytopeCSC format"}},
                                     {"3",{3, "PolytopeJIT format"}},
                                     {"4",{4, "Polyvest: alternative lib, only for single body polytopes - will preprocess first!"}},
                                    }));

    bool doPreprocess = false;
    cliFun.add_long(new CLIF_Option<bool>(&doPreprocess,'p',"doPreprocess","false", {
                                                     {"false",{false,"no preprocessing (may assert)."}},
						     {"true", {true, "Preprocess before running benchmark (could take a while)."}} }));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    
    //Solved_Body* solved_body = solved_body_generator()->get(generator,false);

    Solved_Body *solved_body = generate_read_vinci_polytope(filename);

    if(doPreprocess) {
        std::cout << "Preprocessing as requested...\n";
        Solved_Body* tmp = solved_body->preprocess();
        delete solved_body;
        solved_body = tmp;
        std::cout << "Preprocessing done.\n";
    }

    if(polytopeOptimize) {
        solved_body->optimize();
    }

    switch(polytopeType) {
    case 1: // column major
        solved_body->polytopeTranspose();
        break;
    case 0: // row major
        break;
    case 2: // CSC format
        solved_body->polytopeCSC();
        break;
    case 3: // JIT format
        solved_body->polytopeJIT();
        break;
    }
    if(printBody) {solved_body->print();}
    assert(solved_body->is_normalized);
    FT r0 = 1.0;
    FT r1 = 2*solved_body->n;
    
    ArbitraryExpNum res = volume(solved_body->n, r0, r1, solved_body->bcount, (const void**)solved_body->body, (const Body_T**)solved_body->type);

    std::cout << "Volume is " << res.num << "\n";
    return 0;    
}
