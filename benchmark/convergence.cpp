#include "../src/volume/volume_helper.hpp"


int main(int argc, char *argv[]){

    CLI cli(argc, argv, "convergence");
    CLIFunctionsVolume cliFun(cli);

    cliFun.claimOpt('b', "benchmark configuration");
    std::string body = "cube_r1.0_10";
    // nr of samples taken per configuration
    int nsamples = 1;
    // body to run convergence benchmark on
    auto &body_map = solved_body_generator()->gen_map();
    cliFun.add(new CLIF_Option<std::string>(&body, 'b', "b", body, body_map));
    cliFun.add(new CLIF_OptionNumber<int>(&nsamples, 'b', "s", "1", 1, 100));


    cliFun.claimOpt('w', "walk_size range");
    // lower and upper bound and step size of walk_size parameters for volume algo
    // -> we try all configurations in range(lower, upper+1, step)
    int walk_lower = 1, walk_upper = 10, walk_step = 2;
    cliFun.add(new CLIF_OptionNumber<int>(&walk_lower, 'w', "a", "1", 1, 1000000));
    cliFun.add(new CLIF_OptionNumber<int>(&walk_upper, 'w', "b", "10", 1, 1000000));
    cliFun.add(new CLIF_OptionNumber<int>(&walk_step, 'w', "c", "2", 2, 1000000));


    cliFun.claimOpt('s', "step_size range");
    // lower and upper bounds and step size of step_size parameter for volume algo
    // -> we try all configurations in range(lower, upper+1, step)
    int step_lower = 1600, step_upper = 10000, step_step = 10;
    cliFun.add(new CLIF_OptionNumber<int>(&step_lower, 's', "a", "1600", 1, 1000000));
    cliFun.add(new CLIF_OptionNumber<int>(&step_upper, 's', "b", "10000", 1, 1000000));
    cliFun.add(new CLIF_OptionNumber<int>(&step_step, 's', "c", "10", 2, 1000000));

    cliFun.preParse();
    if (!cli.parse()) {return -1;}
    cliFun.postParse();

    std::cout << "Getting body: " << body << std::endl;
    Solved_Body *b = solved_body_generator()->get(body);

    if (b->volume) {
        std::cout << "Solved volume is: " << b->volume << std::endl;
    }
    else {
        std::cout << "Solved volume is unknown, cannot procede with analysis!\n";
        return 1;
    }
    
    if (!b->is_normalized){
        std::cout << "Preprocess body as it isn't normalized yet...";
        b = b->preprocess();
        std::cout << " done!" << std::endl;
    }
    
    for (int s = step_lower; s <= step_upper; s *= step_step) {
        std::cout << ",\t" << s;
    }

    
    // loop over range of values for step_size and walk_size
    for (int w = walk_lower; w <= walk_upper; w *= walk_step){
        std::cout << std::endl << w;
        for (int s = step_lower; s <= step_upper; s *= step_step){
            step_size = s;
            walk_size = w;

            FT rel_eps = 0;
            for (int i = 0; i < nsamples; i++){
                ArbitraryExpNum vol = volume(b->n, 1.0, (FT) 2*b->n, b->bcount,(const void **) b->body, (const Body_T **) b->type);

                rel_eps += abs(vol.num - b->volume)/b->volume;
            }

            std::cout << ",\t" << rel_eps/nsamples;
        }
    }

    std::cout << std::endl;             
       

    
}
