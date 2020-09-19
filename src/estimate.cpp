#include "volume/volume_helper.hpp"



int main(int argc, char *argv[]){
    CLI_LONG cli(argc,argv,"main");
    CLI_LONG_Functions cliFun(cli);
    initVolumeFunctions(cliFun);
    
    std::string filename = "";
    cliFun.add_long(new CLIF_MandatoryString(&filename, 'f', "filename", "Input filename, file of .ine format."));

    
    bool polytopeOptimize = false; // does not do much anyway, so just drop it
    //cliFun.add_long(new CLIF_Option<bool>(&polytopeOptimize,'o',"polytopeOptimize","false", {
    //                                                 {"false",{false,"-"}},
    // 						     {"true", {true, "Sort constraints to optimize access pattern"}} }));

    bool printBody = false;
    cliFun.add_long(new CLIF_Option<bool>(&printBody,0,"printBody","false", {
                                                     {"false",{false,"-"}},
						     {"true", {true, "Print body."}} },
						     "If true: body is printed to stdout after preprocessing and before estimation is performed."));
    
    double densityThreashold = 0.1;
    cliFun.add_long(new CLIF_OptionNumber<double>(&densityThreashold,'d',"densityThreashold","0.1", 0, 1,
			    "Compared to density after preprocessing. If higher: dense, if lower: sparse."));

    int sparseType = 0;
    cliFun.add_long(new CLIF_Option<int>(&sparseType,0,"sparseType","CSC",
                                    {
                                     {"CSC",{2, "PolytopeCSC format"}},
                                     {"JIT",{3, "PolytopeJIT format"}},
				     },
				     "Data type for sparse case estimation."));

    int denseType = 1;

    bool doPreprocess = false;
    cliFun.add_long(new CLIF_Option<bool>(&doPreprocess,'p',"doPreprocess","false", {
                                                     {"false",{false,"no preprocessing (may assert)."}},
						     {"true", {true, "Preprocess before running benchmark (could take a while)."}} },
						     "TODO: make this isNormalized, a promise. if false -> perform preprocessing"));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();
    
    prng_init();
    
    // we now only load from file
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

    // check density:
    double density = solved_body->density();
    std::cout << "Body Density after preprocessing: " << density << " (vs threashold: "<< densityThreashold<<")\n";
    int polytopeType = 1;
    if(density <= densityThreashold) {
       // sparse
       polytopeType = sparseType;
    } else {
       // dense
       polytopeType = denseType;
    }

    switch(polytopeType) {
    case 1: // column major
	std::cout << "Estimation on PolytopeT Datatype.\n";
        solved_body->polytopeTranspose();
        break;
    case 0: // row major
	std::cout << "Estimation on Polytope Datatype.\n";
        break;
    case 2: // CSC format
	std::cout << "Estimation on PolytopeCSC sparse Datatype.\n";
        solved_body->polytopeCSC();
        break;
    case 3: // JIT format
	std::cout << "Estimation on PolytopeJIT sparse Datatype.\n";
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
